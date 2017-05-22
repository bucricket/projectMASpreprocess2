#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:38:31 2017

@author: mschull
"""
import os
import subprocess
import tarfile
import numpy as np
import glob
from osgeo import gdal,osr
import pandas as pd
from numba import jit
import urllib2, base64

def folders(base):
    dataBase = os.path.join(base,'data')
    landsatDataBase = os.path.join(dataBase,'Landsat-8')
    metBase = os.path.join(dataBase,'MET')
    if not os.path.exists(metBase):
        os.makedirs(metBase) 
    ALEXIbase = os.path.join(dataBase,'ALEXI')
    if not os.path.exists(ALEXIbase):
        os.makedirs(ALEXIbase) 
    resultsBase = os.path.join(base,'outputs')
    albedoBase = os.path.join(landsatDataBase,'albedo')
    if not os.path.exists(albedoBase):
        os.makedirs(albedoBase)   
    landsatSR = os.path.join(landsatDataBase,'SR')
    if not os.path.exists(landsatSR):
        os.makedirs(landsatSR)
    if not os.path.exists(resultsBase):
        os.makedirs(resultsBase)
    landsatLC = os.path.join(landsatDataBase,'LC')
    if not os.path.exists(landsatLC):
        os.makedirs(landsatLC)
    out = {'dataBase':dataBase,'metBase':metBase,'ALEXIbase':ALEXIbase,
           'landsatDataBase':landsatDataBase,'resultsBase':resultsBase,
           'landsatLC':landsatLC,'albedoBase':albedoBase,'landsatSR':landsatSR}
    return out
    
def warp(args):
    """with a def you can easily change your subprocess call"""
    # command construction with binary and options
    options = ['gdalwarp']
    options.extend(args)
    # call gdalwarp 
    subprocess.check_call(options)

def writeArray2Tiff(data,res,UL,inProjection,outfile,outFormat):

    xres = res[0]
    yres = res[1]

    ysize = data.shape[0]
    xsize = data.shape[1]

    ulx = UL[0] #- (xres / 2.)
    uly = UL[1]# - (yres / 2.)
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, xsize, ysize, 1, outFormat)
    #ds = driver.Create(outfile, xsize, ysize, 1, gdal.GDT_Int16)
    
    srs = osr.SpatialReference()
    
    if isinstance(inProjection, basestring):        
        srs.ImportFromProj4(inProjection)
    else:
        srs.ImportFromEPSG(inProjection)
        
    ds.SetProjection(srs.ExportToWkt())
    
    gt = [ulx, xres, 0, uly, 0, -yres ]
    ds.SetGeoTransform(gt)
    
    ds.GetRasterBand(1).WriteArray(data)
    #ds = None
    ds.FlushCache()    
    
    
def getParFromExcel(data,landsatLC,classification,varName):
    ''' Maps LC classification based variables

    Parameters
    ----------
    data : int32
        classification map
    classification : string
        name of LC classification scheme excel tab (i.e. NLCD)
    varName : string
        name of variable

    Returns
    -------
    outVarArray : float
        Mapped vaiable based on LC classification
    '''
    lc = pd.ExcelFile(os.path.join(landsatLC,'landcover.xlsx'))
    lcDF = lc.parse(classification)   
    LCdata = data
    if data.ndim==1:
        outVarArray = np.zeros((data.shape[0]), dtype=np.float)
    else:
        outVarArray = np.zeros((data.shape[0],data.shape[1]), dtype=np.float)
    for row in lcDF.itertuples():
        if classification=='NLCD':
            outVarArray[LCdata == eval('row.%s' % 'NLCD_class')]=eval('row.%s' % varName)
        else:
            outVarArray[LCdata == eval('row.%s' % 'Class')]=eval('row.%s' % varName)
    return outVarArray

def convertBin2tif(inFile,inUL,shape,res):
    inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    outFormat = gdal.GDT_Float32
    read_data = np.fromfile(inFile, dtype=np.float32)
    dataset = np.flipud(read_data.reshape([shape[0],shape[1]]))
    outTif = inFile[:-4]+".tif"
    writeArray2Tiff(dataset,res,inUL,inProj4,outTif,outFormat) 
    
def km2deg(x,y,lat):
    
    # how many KM in 1 deg
    degLat = 110.54 # KM
    degLon = 111.320*np.cos(np.deg2rad(lat))  #KM
    
    degOut = []
    degOut.append(y/degLat)
    degOut.append(x/degLon)
    
    return degOut 
    
def untar(fname, fpath):
    if (fname.endswith('tar.gz') or fname.endswith('tar.bz')):
        tar = tarfile.open(fname)
        tar.extractall(path = fpath)
        tar.close()
        os.remove(fname)
        
def buildvrt(cmd):
    import shlex
    """with a def you can easily change your subprocess call"""
    args = shlex.split(cmd)
    args = args[:-1] + glob.glob(args[-1])
    # This should work now
    subprocess.call(args)

def translate(cmd):
    import shlex
    """with a def you can easily change your subprocess call"""
    args = shlex.split(cmd)
    p = subprocess.call(args)

def clean(directory,fileString):
    from path import path
    d = path(directory)
    files = d.walkfiles(fileString)
    for file in files:
        file.remove()
        print "Removed {} file".format(file)
 
@jit(['float64[:,:](float64[:,:],float64,float64)'])  
def interpOverpassHour(dataset,overpassTime,hours=24.):
    numPixs = dataset.shape[1]
    stack = np.empty([1,numPixs])
    stack[:]=np.nan
    
    #stackReshp = np.reshape(stack,[24,600*1440])
    for j in xrange(numPixs):
        y = dataset[:,j]
        if np.sum(y)==0:
            stack[:,j]=0.0
        else:
            x = range(0,int(hours),int(hours/dataset.shape[0]))
            newX = xrange(int(hours))
            newX = overpassTime
            
            #f = interp1d(x,y, kind='cubic')
            stack[:,j]=np.interp(newX,x,y)
        #stack[:,i]=f(newX)
    
    return stack

# helper function
def _test_outside(testx, lower, upper):
    """
    True if testx, or any element of it is outside [lower, upper].

    Both lower bound and upper bound included
    Input: Integer or floating point scalar or Numpy array.
    """
    test = np.array(testx)
    return np.any(test < lower) or np.any(test > upper)

# custom exception
class RasterError(Exception):
    """Custom exception for errors during raster processing in Pygaarst"""
    pass

class earthDataHTTPRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return urllib2.HTTPRedirectHandler.http_error_302(self, req, fp, code, msg, headers)
    

def getHTTPdata(url,outFN,auth=None):
    request = urllib2.Request(url) 
    if not (auth == None):
        username = auth[0]
        password = auth[1]
        base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')
        request.add_header("Authorization", "Basic %s" % base64string) 
    
    cookieprocessor = urllib2.HTTPCookieProcessor()
    opener = urllib2.build_opener(earthDataHTTPRedirectHandler, cookieprocessor)
    urllib2.install_opener(opener) 
    r = opener.open(request)
    result = r.read()
    
    with open(outFN, 'wb') as f:
        f.write(result)
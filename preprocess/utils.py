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

def folders(base):
    inputDataBase = os.path.join(os.sep,'data','data123','chain','GETD_FINAL')
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
    ndviBase = os.path.join(landsatDataBase,'ndvi')
    if not os.path.exists(ndviBase):
        os.makedirs(ndviBase)
    landsatSR = os.path.join(landsatDataBase,'SR')
    if not os.path.exists(landsatSR):
        os.makedirs(landsatSR)
    if not os.path.exists(resultsBase):
        os.makedirs(resultsBase)
    landsatDN = os.path.join(landsatDataBase,'DN')
    if not os.path.exists(landsatDN):
        os.makedirs(landsatDN)
    landsatLC = os.path.join(landsatDataBase,'LC')
    if not os.path.exists(landsatLC):
        os.makedirs(landsatLC)
    out = {'dataBase':dataBase,'metBase':metBase,'inputDataBase':inputDataBase,
    'landsatDN':landsatDN,'ALEXIbase':ALEXIbase,'landsatDataBase':landsatDataBase,
    'resultsBase':resultsBase,'landsatLC':landsatLC,'albedoBase':albedoBase,
    'ndviBase':ndviBase,'landsatSR':landsatSR}
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
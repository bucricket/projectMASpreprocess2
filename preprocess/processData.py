#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:30:05 2017

@author: mitchell schull
"""
from .landsatTools import landsat_metadata,GeoTIFF
import os
from datetime import datetime
import subprocess
import numpy as np
import utm
import shutil
from .utils import writeArray2Tiff,warp,folders
from .utils import getHTTPdata
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly
from pydap.client import open_url
from pydap.cas import urs
import pygrib
import zipfile
from xml.dom import minidom
from xml.etree.ElementTree import parse, Element, SubElement, tostring
import glob
import gzip


def tile2latlon(tile):
    row = tile/24
    col = tile-(row*24)
    # find lower left corner
    lat= (75.-row*15.)-15.
    lon=(col*15.-180.)-15. 
    return [lat,lon]

def latlon2tile(lat,lon):
    row = int((75-lat)/15)
    col = int((abs(-180-lon)/15)+1)  
    tile = (row*24+col)
    return tile
   
class Landsat(object):
    def __init__(self, filepath,inputLC):
#        base = os.path.abspath(os.path.join(filepath,os.pardir,os.pardir,os.pardir,
#                                            os.pardir,os.pardir))
        base = os.getcwd()
        Folders = folders(base)    
        self.landsatLC = Folders['landsatLC']
        self.satCache = Folders['satCache']
        self.inputLC = inputLC
        meta = landsat_metadata(filepath)
        self.sceneID = meta.LANDSAT_SCENE_ID
        self.productID = filepath.split(os.sep)[-1][:-8]
        self.scene = self.sceneID[3:9]
        if meta.SPACECRAFT_ID == 'LANDSAT_7':        
            self.sat = 7
        else:
            self.sat = 8  
        satscene_path = os.path.join(self.satCache,"LANDSAT","L%d" % self.sat,self.scene)
        self.lc_path = os.path.join(satscene_path,"LC")
        if not os.path.exists(self.lc_path):
            os.mkdir(self.lc_path)
        self.albedo_path = os.path.join(satscene_path,"ALBEDO")
        if not os.path.exists(self.albedo_path):
            os.mkdir(self.albedo_path)
        self.landsatSR = os.path.join(satscene_path,"RAW_DATA")
#        self.productID = meta.LANDSAT_PRODUCT_ID
        ls = GeoTIFF(os.path.join(self.landsatSR,'%s_sr_band1.tif' % self.productID))
        self.nrow = ls.nrow
        self.ncol = ls.ncol
        self.proj4 = ls.proj4
        self.inProj4 = ls.proj4
        self.ulx = meta.CORNER_UL_PROJECTION_X_PRODUCT
        self.uly = meta.CORNER_UL_PROJECTION_Y_PRODUCT
        self.lrx = meta.CORNER_LR_PROJECTION_X_PRODUCT
        self.lry = meta.CORNER_LR_PROJECTION_Y_PRODUCT
        self.ulLat = meta.CORNER_UL_LAT_PRODUCT
        self.ulLon = meta.CORNER_UL_LON_PRODUCT
        self.lrLat = meta.CORNER_LR_LAT_PRODUCT
        self.lrLon = meta.CORNER_LR_LON_PRODUCT
        self.delx = meta.GRID_CELL_SIZE_REFLECTIVE
        self.dely = meta.GRID_CELL_SIZE_REFLECTIVE

        
    def getLC(self,classification):
#        scene = self.scene
        sceneID = self.sceneID
        if classification=='NLCD': #NLCD
            outfile = os.path.join(self.landsatLC,'nlcd_2011_landcover_2011_edition_2014_10_10','nlcd_2011_landcover_2011_edition_2014_10_10.img')
            if not os.path.exists(outfile):
                print "download NLCD data at:\n"
                print 'http://www.landfire.gov/bulk/downloadfile.php?TYPE=nlcd2011&FNAME=nlcd_2011_landcover_2011_edition_2014_10_10.zip'
                print "and unzip in %s\n" % self.landsatLC
                raise SystemExit(0)
            self.inProj4 = '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
        else:
            LCtemp = os.path.join(self.landsatLC,"temp")
            if not os.path.exists(LCtemp):
                os.mkdir(LCtemp)                                 
            #get UTM info
            utmZone= []
            corners = [[self.ulLat,self.ulLon],[self.ulLat,self.lrLon],
                       [self.lrLat,self.lrLon],[self.lrLat,self.ulLon]]
            for i in xrange(4):
                lat =corners[i][0]
                lon =corners[i][1]
                utmZone.append(utm.from_latlon(lat,lon)[2])
            latNames = [np.floor(self.ulLat/5.)*5,np.floor(self.lrLat/5.)*5.]
            utmzone = np.unique(utmZone)
            latNames = np.unique(latNames)
            for i in xrange(len(utmzone)):
                for j in xrange(len(latNames)):
                    if latNames[j]>=0:
                        hemisphere = 'N'
                    else:
                        hemisphere = 'S'
                    LCdataFolder  = os.path.join(self.inputLC,'%s%d_%d_2010LC030' % (hemisphere,utmzone[i],latNames[j]))
                    if not os.path.exists(LCdataFolder):
                        zipFN = os.path.join("%s.zip" % LCdataFolder)
                        zip_ref = zipfile.ZipFile(zipFN, 'r')
                        zip_ref.extractall(self.inputLC)
                        zip_ref.close()
                    LCdata = os.path.join(LCdataFolder,'%s%d_%d_2010lc030.tif' % (hemisphere.lower(),utmZone[i],latNames[j]))
                    inputLCdata = os.path.join(LCtemp,LCdata.split(os.sep)[-1])
                    if not os.path.exists(inputLCdata):
                        os.symlink(LCdata,inputLCdata)

            # mosaic dataset if needed
            outfile = os.path.join(LCtemp,'tempMos.tif')

            subprocess.check_output('gdalbuildvrt -srcnodata 0 %s.vrt %s%s*.tif' % (outfile[:-4], LCtemp,os.sep),shell=True)
            subprocess.call(["gdal_translate", "-of", "GTiff", "%s.vrt" % outfile[:-4],"%s" % outfile])

#            #====remove unzipped folders
#            LCfolders=next(os.walk(self.inputLC))[1]
#            for LCfolder in LCfolders:
#                shutil.rmtree(os.path.join(self.inputLC,LCfolder))

        outfile2=os.path.join(self.lc_path,'%s%s' % (sceneID,'_LC.tiff'))

        optionList = ['-overwrite', '-s_srs', '%s' % self.inProj4,'-t_srs','%s' % self.proj4,\
        '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,\
        '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % outfile, '%s' % outfile2]
        warp(optionList)
        shutil.rmtree(LCtemp)
        


    def getAlbedo(self):
#        scene = self.scene
        sceneID = self.sceneID
#        sceneFolderAlbedo = os.path.join(self.albedoBase,scene)
#        if not os.path.exists(sceneFolderAlbedo):
#            os.makedirs(sceneFolderAlbedo)
#        albedoPath = os.path.join(sceneFolderAlbedo,'%s_albedo.tiff' % sceneID)
        albedoPath = os.path.join(self.albedo_path,'%s_albedo.tiff' % sceneID)
        
    
        bands = [1,3,4,5,7] # dont use blue
        # extract the desired surface refelctance bands from landsat
        data = []   
        for i in xrange(len(bands)):
            landsatgeotif = os.path.join(self.landsatSR,'%s_sr_band%d.tif' % (self.productID,bands[i]))
            if os.path.exists(landsatgeotif):
                ls = GeoTIFF(landsatgeotif)
                data.append(ls.data*0.0001)
            else:
                print "no sr file for scenedID: %s, band: %d" % (sceneID,bands[i])
        albedotemp=(0.356*data[0])+(0.130*data[1])+(0.373*data[2])+(0.085*data[3])+(0.072*data[4])-0.0018
        ls.clone(albedoPath,albedotemp)
        
class ALEXI:
    def __init__(self, filepath):
#        base = os.path.abspath(os.path.join(filepath,os.pardir,os.pardir,os.pardir,
#                                            os.pardir,os.pardir))
        base = os.getcwd()
        Folders = folders(base)    
        self.ALEXIbase = Folders['ALEXIbase']
        self.satCache = Folders['satCache']
        self.modCache = Folders['modCache']
        meta = landsat_metadata(filepath)
        self.sceneID = meta.LANDSAT_SCENE_ID
#        self.productID = meta.LANDSAT_PRODUCT_ID
        self.productID = filepath.split(os.sep)[-1][:-8]
        self.scene = self.sceneID[3:9]
        if meta.SPACECRAFT_ID == 'LANDSAT_7':        
            self.sat = 7
        else:
            self.sat = 8  
        satscene_path = os.path.join(self.satCache,"LANDSAT","L%d" % self.sat,self.scene)
        self.alexi_path = os.path.join(satscene_path,"ET","400m")
        if not os.path.exists(self.alexi_path):
            os.mkdir(self.alexi_path)
        self.landsatSR = os.path.join(satscene_path,"RAW_DATA")
        self.yeardoy = self.sceneID[9:16]
        ls = GeoTIFF(os.path.join(self.landsatSR,'%s_sr_band1.tif' % self.productID))
        self.proj4 = ls.proj4
        self.nrow = ls.nrow
        self.ncol = ls.ncol
        self.ulx = meta.CORNER_UL_PROJECTION_X_PRODUCT
        self.uly = meta.CORNER_UL_PROJECTION_Y_PRODUCT
        self.lrx = meta.CORNER_LR_PROJECTION_X_PRODUCT
        self.lry = meta.CORNER_LR_PROJECTION_Y_PRODUCT
        self.ulLat = meta.CORNER_UL_LAT_PRODUCT
        self.ulLon = meta.CORNER_UL_LON_PRODUCT
        self.lrLat = meta.CORNER_LR_LAT_PRODUCT
        self.lrLon = meta.CORNER_LR_LON_PRODUCT
        self.delx = meta.GRID_CELL_SIZE_REFLECTIVE
        self.dely = meta.GRID_CELL_SIZE_REFLECTIVE
        self.landsatDate = meta.DATE_ACQUIRED
        self.landsatTime = meta.SCENE_CENTER_TIME[:-2]
        d = datetime.strptime('%s%s' % (self.landsatDate,self.landsatTime),'%Y-%m-%d%H:%M:%S.%f')
        self.hr = d.hour #UTC
        self.year = d.year
        self.inputET = os.path.join(self.modCache,'ALEXI')
#        self.inputET = inputET
        
        
    def getALEXIdata(self,ALEXIgeodict,isUSA):       
        #ALEXI input info
#        ALEXI_ulLon = ALEXIgeodict['ALEXI_ulLon']
#        ALEXI_ulLat = ALEXIgeodict['ALEXI_ulLat']
        inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        ETtemp = os.path.join(self.ALEXIbase,"temp")
        if not os.path.exists(ETtemp):
            os.makedirs(ETtemp)
            
        outfile = os.path.join(self.alexi_path,'%s_alexiET.tiff' % self.sceneID)
        if not os.path.exists(outfile):
            print 'processing : %s...' % outfile
#            subsetFile = outfile[:-5]+'Sub.tiff'
            #********THIS SECTION IS A TEMPERARY FIX
            corners = [[self.ulLat,self.ulLon],[self.ulLat,self.lrLon],
                       [self.lrLat,self.lrLon],[self.lrLat,self.ulLon]]
            tile_num =[]
#            ULlat =[]
#            ULlon =[]
            for i in xrange(4):
                lat =corners[i][0]
                lon =corners[i][1]
                tile  = latlon2tile(lat,lon)
#                LLlat,LLlon = tile2latlon(tile)
#                row = int((75-lat)/15)
#                col = int((abs(-180-lon)/15)+1)
#                ULlat.append(LLlat+15.0)
#                ULlon.append(LLlon)      
                tile_num.append(tile)

            for i in xrange(len(tile_num)):
#                inUL = [ULlon[i],ULlat[i]]
#                ETdata = os.path.join(self.inputET,
#                                      'FINAL_EDAY_%s_T%03d.dat' % (int(self.sceneID[9:16]),tile_num[i]))
                ETdata = os.path.join(self.inputET,
                                      'FINAL_EDAY_%s_T%03d.tif' % (int(self.sceneID[9:16]),tile_num[i]))
                localETpath = os.path.join(ETtemp,ETdata.split(os.sep)[-1])
                if not os.path.exists(os.path.join(ETtemp,localETpath)):
                    if not os.path.exists(ETdata):
                        print("data doesn't exist!")
                        continue
                    else:
                        os.symlink(ETdata,os.path.join(ETtemp,localETpath))
#                convertBin2tif(localETpath,inUL,ALEXIshape,inRes)
#                os.remove(os.path.join(ETtemp,localETpath))
             # mosaic dataset if needed
            outfile2 = os.path.join(self.ALEXIbase,'tempMos.tif')

            subprocess.check_output('gdalbuildvrt %s.vrt %s%s*.tif' % (outfile2[:-4], ETtemp,os.sep),shell=True)
            subprocess.call(["gdal_translate", "-of", "GTiff", "%s.vrt" % outfile2[:-4],"%s" % outfile2])

            #========fill in missing data from VIIRS and Landsat data======
            g = gdal.Open(outfile2,GA_ReadOnly)
            et= g.ReadAsArray()
            et[et==-9999]=0
            ls = GeoTIFF(outfile2)
            mask = os.path.join(ETtemp,"Mask.tif")
            masked = os.path.join(ETtemp,"Masked.tif")
            ls.clone(mask,et)
            subprocess.check_output('gdal_fillnodata.py %s %s -mask %s -of GTiff' % (outfile2,masked,mask),shell=True)
            
            optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
            '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'near',\
            '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % masked , '%s' % outfile]
            warp(optionList)
            shutil.rmtree(ETtemp)
class MET:
    def __init__(self, filepath,session):
#        base = os.path.abspath(os.path.join(filepath,os.pardir,os.pardir,os.pardir,
#                                            os.pardir,os.pardir))
        base = os.getcwd() 
        Folders = folders(base)
        self.earthLoginUser = session[0]
        self.earthLoginPass = session[1]
        self.metBase = Folders['metBase']
        self.satCache = Folders['satCache']
        meta = landsat_metadata(filepath)
        self.sceneID = meta.LANDSAT_SCENE_ID
#        self.productID = meta.LANDSAT_PRODUCT_ID
        self.productID = filepath.split(os.sep)[-1][:-8]
        self.scene = self.sceneID[3:9]
        if meta.SPACECRAFT_ID == 'LANDSAT_7':        
            self.sat = 7
        else:
            self.sat = 8  
        satscene_path = os.path.join(self.satCache,"LANDSAT","L%d" % self.sat,self.scene)
        self.gsip_path = os.path.join(self.satCache,"GSIP")
        if not os.path.exists(self.gsip_path):
            os.mkdir(self.gsip_path)
        self.met_path = os.path.join(satscene_path,"MET")
        if not os.path.exists(self.met_path):
            os.mkdir(self.met_path)
        self.insol_path = os.path.join(satscene_path,"INSOL")
        if not os.path.exists(self.insol_path):
            os.mkdir(self.insol_path)
        self.landsatSR = os.path.join(satscene_path,"RAW_DATA")
        self.yeardoy = self.sceneID[9:16]
        ls = GeoTIFF(os.path.join(self.landsatSR,'%s_sr_band1.tif' % self.productID))
        self.proj4 = ls.proj4
        self.nrow = ls.nrow
        self.ncol = ls.ncol
        self.ulx = meta.CORNER_UL_PROJECTION_X_PRODUCT
        self.uly = meta.CORNER_UL_PROJECTION_Y_PRODUCT
        self.lrx = meta.CORNER_LR_PROJECTION_X_PRODUCT
        self.lry = meta.CORNER_LR_PROJECTION_Y_PRODUCT
        self.ulLat = meta.CORNER_UL_LAT_PRODUCT
        self.ulLon = meta.CORNER_UL_LON_PRODUCT
        self.lrLat = meta.CORNER_LR_LAT_PRODUCT
        self.lrLon = meta.CORNER_LR_LON_PRODUCT
        self.delx = meta.GRID_CELL_SIZE_REFLECTIVE
        self.dely = meta.GRID_CELL_SIZE_REFLECTIVE
        self.landsatDate = meta.DATE_ACQUIRED
        self.landsatTime = meta.SCENE_CENTER_TIME[:-2]
        self.d = datetime.strptime('%s%s' % (self.landsatDate,self.landsatTime),'%Y-%m-%d%H:%M:%S.%f')
        d = self.d
        self.year = d.year
        self.month = d.month
        self.day = d.day
        self.hr = d.hour #UTC
        
                
    def getCFSR(self):
        
        # set Met dataset geo info
        CFSR_ulLat = 90.0
        CFSR_ulLon = -180.0
        CFSRLatRes = 0.50
        CFSRLonRes = 0.50
        inRes = [CFSRLonRes,CFSRLatRes]
        inUL = [CFSR_ulLon,CFSR_ulLat]
        inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        
        #===Open CFSR file=================
        dailyPath = os.path.join(self.metBase,'%s' % self.scene) 
        if not os.path.exists(dailyPath):
            os.makedirs(dailyPath)
        ncdcURL = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_pgbh/'
        dateDiff = (datetime.today()-self.d).days
        if (dateDiff<7):
            ncdcURL = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/cfs/prod/'

        iHour = (int(self.hr/6.)*6.)
        fHour = self.hr-iHour
        hr1file = 'cdas1.t%02dz.pgrbh%02d.grib2' % (iHour,fHour)
        pydapURL = os.path.join(ncdcURL,"%s" % self.year,"%d%02d" % (self.year,self.month),"%d%02d%02d" % (self.year,self.month,self.day),hr1file)
        outFN = os.path.join(os.getcwd(),hr1file)
        getHTTPdata(pydapURL,outFN)
        print pydapURL
        
        #=====Surface Pressure=============        
        grbs = pygrib.open(outFN)

        sfcP = []
        for grb in grbs:
            if grb.name == 'Surface pressure' and grb.typeOfLevel == 'surface': 
                sfcP.append(grb.values)
        opdata = sfcP[0]
        #---flip data so ul -180,90-----------      
        b = opdata[:,:360]
        c = opdata[:,360:]
        flipped = np.hstack((c,b))
        print np.shape(flipped)
        #----write out global dataset to gTiff--------
        outfile = os.path.join(dailyPath,'%s_p.tiff' % (self.sceneID))
        outFormat = gdal.GDT_Float32
        writeArray2Tiff(flipped.astype(float),inRes,inUL,inProj4,outfile,outFormat)
        
#        subsetFile = outfile[:-5]+'Sub.tiff'
        subsetFile = os.path.join(self.met_path,'%s_p.tiff' % (self.sceneID))
        optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
                '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'bilinear',\
                '-ts', '%f' % self.ncol, '%f' % self.nrow ,'-multi','-of','GTiff','%s' % outfile, '%s' % subsetFile]
        print "running pydisalexi from Aug. 11"
        warp(optionList)   
        
        #====Wind Speed================        
        grbs.rewind()
        Ucomp = []
        for grb in grbs:
            if grb.name == 'U component of wind' and grb.typeOfLevel == 'sigma': 
                Ucomp.append(grb.values)
        Uopdata = Ucomp[0]
        
        grbs.rewind()
        Vcomp = []
        for grb in grbs:
            if grb.name == 'V component of wind' and grb.typeOfLevel == 'sigma': 
                Vcomp.append(grb.values)
        Vopdata = Vcomp[0]
        
        wind = np.sqrt(Uopdata**2+Vopdata**2)
        #---flip data so ul -180,90-----------      
        b = wind[:,:360]
        c = wind[:,360:]
        flipped = np.hstack((c,b))
        #----write out global dataset to gTiff--------
        
        outfile = os.path.join(dailyPath,'%s_u.tiff' % (self.sceneID))
        outFormat = gdal.GDT_Float32
        writeArray2Tiff(flipped.astype(float),inRes,inUL,inProj4,outfile,outFormat)
        
#        subsetFile = outfile[:-5]+'Sub.tiff'
        subsetFile = os.path.join(self.met_path,'%s_u.tiff' % (self.sceneID))
        optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
                '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'bilinear',\
                '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % outfile, '%s' % subsetFile]
        
        warp(optionList)
        
        #===== Specific humidity =============
        grbs.rewind()
        specH = []
        for grb in grbs:
            if grb.name == 'Specific humidity' and grb.typeOfLevel == 'hybrid': 
                specH .append(grb.values)
        opdata = specH [0]
        #---flip data so ul -180,90-----------      
        b = opdata[:,:360]
        c = opdata[:,360:]
        flipped = np.hstack((c,b))
        #----write out global dataset to gTiff--------
        
        outfile = os.path.join(dailyPath,'%s_q2.tiff' % (self.sceneID))
        outFormat = gdal.GDT_Float32
        writeArray2Tiff(flipped.astype(float),inRes,inUL,inProj4,outfile,outFormat)
        
#        subsetFile = outfile[:-5]+'Sub.tiff'
        subsetFile = os.path.join(self.met_path,'%s_q2.tiff' % (self.sceneID))
        optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
                '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'bilinear',\
                '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % outfile, '%s' % subsetFile]
        warp(optionList)     
        
        #===== Temperature 2m =============
        grbs.rewind()
        Temp = []
        for grb in grbs:
            if grb.name == 'Temperature' and grb.typeOfLevel == 'sigma': 
                Temp.append(grb.values)
        opdata = Temp[0]
        #---flip data so ul -180,90-----------      
        b = opdata[:,:360]
        c = opdata[:,360:]
        flipped = np.hstack((c,b))
        #----write out global dataset to gTiff--------
        
        outfile = os.path.join(dailyPath,'%s_Ta.tiff' % (self.sceneID))
        outFormat = gdal.GDT_Float32
        writeArray2Tiff(flipped.astype(float),inRes,inUL,inProj4,outfile,outFormat)
        
#        subsetFile = outfile[:-5]+'Sub.tiff'
        subsetFile = os.path.join(self.met_path,'%s_Ta.tiff' % (self.sceneID))
        optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
                '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'bilinear',\
                '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % outfile, '%s' % subsetFile]
        warp(optionList) 
        shutil.rmtree(dailyPath)
                
    def getInsolation(self):
        MERRA2_ulLat = 90.0
        MERRA2_ulLon = -180.0
        MERRA2LatRes = 0.5
        MERRA2LonRes = 0.625
        inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        #========get MERRA2 Insolation data at overpass time====================== 
        inRes = [MERRA2LonRes,MERRA2LatRes]
        inUL = [MERRA2_ulLon,MERRA2_ulLat]
    
        if self.year <1992:
            fileType = 100
        elif self.year >1991 and self.year < 2001:
            fileType=200
        elif self.year > 2000 and self.year<2011:
            fileType = 300
        else:
            fileType = 400
    
        opendap_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/'
        product = 'M2T1NXRAD.5.12.4'
        
        filename = 'MERRA2_%d.tavg1_2d_rad_Nx.%04d%02d%02d.nc4' % (fileType,self.year,self.month,self.day)
        fullUrl =os.path.join(opendap_url,product,'%04d'% self.year,'%02d'% self.month,filename)
    
        dailyPath = os.path.join(self.metBase,'%s' % self.scene)
        if not os.path.exists(dailyPath):
            os.makedirs(dailyPath)
        #====get overpass hour insolation=========================================
        outFN = os.path.join(self.insol_path,'%s_Insol1Sub.tiff' % self.sceneID) 
        print 'processing : %s...' % outFN
        session = urs.setup_session(username = self.earthLoginUser, 
                    password = self.earthLoginPass,
                    check_url=fullUrl)
        d = open_url(fullUrl,session=session)
        Insol = d.SWGDNCLR
        if not os.path.exists(outFN):
        # wv_mmr = 1.e-6 * wv_ppmv_layer * (Rair / Rwater)
        # wv_mmr in kg/kg, Rair = 287.0, Rwater = 461.5
            dataset = np.squeeze(Insol[self.hr,:,:,:])*1.
            
            outfile = os.path.join(dailyPath,'%s_%s.tiff' % (self.sceneID,'Insol1'))
#            subsetFile = outfile[:-5]+'Sub.tiff'
            subsetFile = os.path.join(self.insol_path,'%s_%s.tiff' % (self.sceneID,'Insol1'))
            outFormat = gdal.GDT_Float32
            writeArray2Tiff(dataset,inRes,inUL,inProj4,outfile,outFormat)
            optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
            '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'bilinear',\
            '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % outfile, '%s' % subsetFile]
            warp(optionList)
            
    #====get daily insolation==================================================
        outFN = os.path.join(self.insol_path,'%s_Insol24Sub.tiff' % self.sceneID)
        if not os.path.exists(outFN):
            dataset2 = np.flipud(np.sum(np.squeeze(Insol[:,:,:]),axis=0))
            outfile = os.path.join(dailyPath,'%s_%s.tiff' % (self.sceneID,'Insol24'))
#            subsetFile = outfile[:-5]+'Sub.tiff'
            subsetFile = os.path.join(self.insol_path,'%s_%s.tiff' % (self.sceneID,'Insol24'))
            outFormat = gdal.GDT_Float32
            writeArray2Tiff(dataset2,inRes,inUL,inProj4,outfile,outFormat)
            optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % self.proj4,\
            '-te', '%f' % self.ulx, '%f' % self.lry,'%f' % self.lrx,'%f' % self.uly,'-r', 'bilinear',\
            '-ts', '%f' % self.ncol, '%f' % self.nrow,'-multi','-of','GTiff','%s' % outfile, '%s' % subsetFile]
            warp(optionList)
            shutil.rmtree(dailyPath)
            
    def getGSIP(self):
        doy = (self.d-datetime(self.year-1,12,31)).days
        date = '%d%03d' % (self.year,doy)
    #=======================GSIP daily insolation==============================
        def gunzip(fn, *positional_parameters, **keyword_parameters):
            inF = gzip.GzipFile(fn, 'rb')
            s = inF.read()
            inF.close()
            if ('out_fn' in keyword_parameters):
                outF = file(keyword_parameters['out_fn'], 'wb')
            else:
                outF = file(fn[:-3], 'wb')
                  
            outF.write(s)
            outF.close()
#        def convertGSIP2tiff(year,doy):
#            inProjection = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
#            date = '%d%03d' % (year,doy)
#            gsip_fn = glob.glob(os.path.join(GSIP_path, "%d" % year,'*.gsipL3_global_GDA_%s.nc.gz' % date))[0]
#            if os.path.exists(gsip_fn):
#                gunzip(gsip_fn)
#                tif_fn = gsip_fn[:-5]+'tif' 
#                nc_fn = 'NETCDF:"%s":insolation' % gsip_fn[:-3]
#                if not os.path.exists(tif_fn):
#                    ds = gdal.Open(nc_fn)
#                    aa = ds.GetRasterBand(1).ReadAsArray()*0.042727217           
#                    writeArray2Tiff(aa,[0.05,0.05],[-180.,90.],inProjection,tif_fn,gdal.GDT_Float32)
#                
#        def processGSIPtiles(tile,year,doy):
#            LLlat,LLlon = tile2latlon(tile)
#            URlat = LLlat+15.
#            LRlon = LLlon+15.
#            date = '%d%03d' % (year,doy)
#            insol24_fn = os.path.join(insol_path, 'RS24_%s_T%03d.tif' % (date,tile))
#            tif_fn = glob.glob(os.path.join(GSIP_path, "%d" % year,'*.gsipL3_global_GDA_%s.tif' % date))[0]
#            if os.path.exists(tif_fn): 
#                outds = gdal.Open(tif_fn)
#                outds = gdal.Translate(insol24_fn, outds,options=gdal.TranslateOptions(xRes=0.004,yRes=0.004,
#                                                                                    projWin=[LLlon,URlat,LRlon,LLlat],
#                                                                                    resampleAlg = 'bilinear'))
#                outds = None
                
    #==============GSIP hourly insolation==================================

    
    
        def addGeoloc(data_fn,lat_fn,lon_fn):
            out_fn = data_fn[:-4]+'.vrt'
            tree = parse(data_fn)
            root = tree.getroot()
            a = root[0]
            b = root[1]
            c = Element('metadata')
            c.set('domain','GEOLOCATION')
            mdi1 = SubElement(c,'mdi')
            mdi1.text = "%s" % lon_fn
            mdi1.set('key','X_DATASET')
            mdi2 = SubElement(c,'mdi')
            mdi2.text = "1"
            mdi2.set('key','X_BAND')
            mdi3 = SubElement(c,'mdi')
            mdi3.text = "%s" % lat_fn
            mdi3.set('key','Y_DATASET')
            mdi4 = SubElement(c,'mdi')
            mdi4.text = "1"
            mdi4.set('key','Y_BAND')
            mdi5 = SubElement(c,'mdi')
            mdi5.text = "0"
            mdi5.set('key','PIXEL_OFFSET')
            mdi6 = SubElement(c,'mdi')
            mdi6.text = "0"
            mdi6.set('key','LINE_OFFSET')
            mdi7 = SubElement(c,'mdi')
            mdi7.text = "1"
            mdi7.set('key','PIXEL_STEP')
            mdi7 = SubElement(c,'mdi')
            mdi7.text = "1"
            mdi7.set('key','LINE_STEP')
            
            root = Element('VRTDataset')
            root.set('rasterXSize',"3712")
            root.set('rasterYSize',"3712")
            root.extend((c,a,b))
        #    xmlstr = minidom.parseString(tostring(root)).toprettyxml(indent="  ")
            
            with open(out_fn, "w") as f:
                f.write(prettify(root))
        
        
        def prettify(elem):
            rough_string = tostring(elem, 'utf-8') #xml as ElementTree
            reparsed = minidom.parseString(rough_string) #mini as minidom
            return '\n'.join([line for line in reparsed.toprettyxml(indent=' '*2).split('\n') if line.strip()]) 
#====convert NETCDF data to VRTs===============================================
#        gsip_path = os.path.join(self.gsip_path, "%d" % self.year) 
#        if not os.path.exists(gsip_path):
#            os.mkdir(gsip_path)
        gsip_path = os.path.join(self.gsip_path,'*gsipL2_met10_MSGFD_%s_%02d30.nc.gz' % (date,self.hr))
        print(gsip_path)
        
        try:
            gsip_fn = glob.glob(gsip_path)[0]    
        except ValueError:
            print "You need to download gsipL2_met10_MSGFD data"
            exit
        layers = ["flux_swd_sfc","pixel_latitude","pixel_longitude"]
        if os.path.exists(gsip_fn):
            gunzip(gsip_fn)
            netcdf_fn = gsip_fn[:-3]
            outFN = os.path.join(self.insol_path,'%s_Insol1Sub.tiff' % self.sceneID)
            if not os.path.exists(outFN):
                fns = []
                for layer in layers:
                    raw_fn = 'NETCDF:%s:%s'  % (netcdf_fn,layer)  
                    vrt_fn = os.path.join(os.getcwd(),'%s.vrt'% layer)
                    outds = gdal.Open(raw_fn)
                    outds = gdal.Translate(vrt_fn, outds,options=gdal.TranslateOptions(format="VRT"))
                    outds = None
                    fns.append(vrt_fn)
                addGeoloc(fns[0],fns[1],fns[2])
                vrt_fn = os.path.join(os.getcwd(),'%s.vrt'% layers[0])                
                in_ds = gdal.Open(vrt_fn)
                outds = gdal.Warp(outFN, in_ds,options=gdal.WarpOptions(resampleAlg='bilinear',geoloc=True,
                                                                         dstSRS="EPSG:4326",
                                                                         outputBounds=(self.ulLon,self.lrLat,self.lrLon,self.ulLat),
                                                                         width=self.ncol,
                                                                         height=self.nrow,
                                                                         multithread=True))
                outds = None
            os.remove(netcdf_fn)
        #=======Daily GSIP====================================================      
#        gsip_fn = glob.glob(os.path.join(self.gsip_path,'*.gsipL3_met10_%s_%d30.nc.gz' % (date,self.hr)))[0]
#        gsip_path = os.path.join(self.gsip_path, "%d" % self.year) 
#        if not os.path.exists(gsip_path):
#            os.mkdir(gsip_path)

#        gsip_fn = glob.glob(os.path.join(self.gsip_path,'*gsipL3_global_GDA_%s.nc.gz' % date))[0]
#        layers = ["insolation","latitude","longitude"]
#        if os.path.exists(gsip_fn):
#            gunzip(gsip_fn)
#            netcdf_fn = gsip_fn[:-3]
#            outFN = os.path.join(self.insol_path,'%s_Insol24Sub.tiff' % self.sceneID)
#            if not os.path.exists(outFN):
#                fns = []
#                for layer in layers:
#                    raw_fn = 'NETCDF:%s:%s'  % (netcdf_fn,layer)  
#                    vrt_fn = os.path.join(os.getcwd(),'%s.vrt'% layer)
#                    outds = gdal.Open(raw_fn)
#                    outds = gdal.Translate(vrt_fn, outds,options=gdal.TranslateOptions(format="VRT"))
#                    outds = None
#                    fns.append(vrt_fn)
#                addGeoloc(fns[0],fns[1],fns[2])
#                vrt_fn = os.path.join(os.getcwd(),'%s.vrt'% layers[0])                
#                in_ds = gdal.Open(vrt_fn)
#                outds = gdal.Warp(outFN, in_ds,options=gdal.WarpOptions(resampleAlg='bilinear',geoloc=True,
#                                                                         dstSRS="EPSG:4326",
#                                                                         outputBounds=(self.ulLon,self.lrLat,self.lrLon,self.ulLat),
#                                                                         width=self.ncol,
#                                                                         height=self.nrow,
#                                                                         multithread=True))
#                outds = None
#            os.remove(netcdf_fn)


        inProjection = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        try:
            gsip_fn = glob.glob(os.path.join(self.gsip_path,'*gsipL3_global_GDA_%s.nc.gz' % date))[0]    
        except ValueError:
            print "You need to download gsipL3_global_GDA data"
            exit
        if os.path.exists(gsip_fn):
            gunzip(gsip_fn)
            tif_fn = gsip_fn[:-5]+'tif' 
            nc_fn = 'NETCDF:"%s":insolation' % gsip_fn[:-3]
            if not os.path.exists(tif_fn):
                ds = gdal.Open(nc_fn)
                aa = ds.GetRasterBand(1).ReadAsArray()*0.042727217           
                writeArray2Tiff(aa,[0.05,0.05],[-180.,90.],inProjection,tif_fn,gdal.GDT_Float32)
            outFN = os.path.join(self.insol_path,'%s_Insol24Sub.tiff' % self.sceneID)
            if not os.path.exists(outFN):
                in_ds = gdal.Open(tif_fn)
                outds = gdal.Warp(outFN, in_ds,options=gdal.WarpOptions(resampleAlg='bilinear',
                                                                         dstSRS="EPSG:4326",
                                                                         outputBounds=(self.ulLon,self.lrLat,self.lrLon,self.ulLat),
                                                                         width=self.ncol,
                                                                         height=self.nrow,
                                                                         multithread=True))
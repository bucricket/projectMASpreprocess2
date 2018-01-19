#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:09:46 2017

@author: mschull
"""
import os
import keyring
import getpass
import argparse
import glob
import pycurl
from .utils import folders
from .processData import ALEXI,MET,Landsat
from .landsatTools import landsat_metadata
from processlai import processlai
from processlst import processlst
from getlandsatdata import getlandsatdata


#===set golbal paths===========================================================
base = os.getcwd()
cacheDir = os.path.abspath(os.path.join(base,os.pardir,"SATELLITE_DATA"))
Folders = folders(base)  
ALEXIbase = Folders['ALEXIbase']
metBase = Folders['metBase']
landsatDataBase = Folders['landsatDataBase']
landsatSR = Folders['landsatSR']

def prepare_data(fn,session,isUSA,LCpath,ETpath):

    LCpath = LCpath
    ETpath = ETpath
    meta = landsat_metadata(fn)
    sceneID = meta.LANDSAT_SCENE_ID
    scene = sceneID[3:9]
    isUSA = isUSA
    ALEXI_ulLon = 0.0 
    ALEXI_ulLat = 0.0
    #--------------------
    ALEXILatRes = 0.004
    ALEXILonRes = 0.004
    ALEXIshape = [3750,3750]
    ALEXIgeodict ={'ALEXI_ulLat':ALEXI_ulLat,'ALEXI_ulLon':ALEXI_ulLon,
               'ALEXI_LatRes':ALEXILatRes,'ALEXI_LonRes':ALEXILonRes,
               'ALEXIshape': ALEXIshape}

        
    #=====pick Landcover map===================================================
    if isUSA ==1:
        landcover = 'NLCD'
    else:
        landcover = 'GlobeLand30'

    
    #=====prepare the ETd data=================================================    
    sceneDir = os.path.join(ALEXIbase,'%s' % scene)        
    
    outFN = os.path.join(sceneDir,'%s_alexiETSub.tiff' % sceneID) 
    if not os.path.exists(outFN):
        print 'get->ALEXI ET...'
        if not os.path.exists(sceneDir):
            os.makedirs(sceneDir)
        a = ALEXI(fn,ETpath)
        a.getALEXIdata(ALEXIgeodict,isUSA)
    
    #=====prepare MET data=====================================================    
    sceneDir = os.path.join(metBase,'%s' % scene)
    
    outFN = os.path.join(sceneDir,'%s_pSub.tiff' % sceneID) 
    if not os.path.exists(outFN):
        print 'get->MET data...'
        if not os.path.exists(sceneDir):
            os.makedirs(sceneDir)
        a = MET(fn,session)
        a.getCFSR()
    
    #====prepare insolation====================================================
    outFN = os.path.join(sceneDir,'%s_Insol1Sub.tiff' % sceneID)
    if not os.path.exists(outFN):
        if not os.path.exists(sceneDir):
            os.makedirs(sceneDir)
        a = MET(fn,session)
        a.getInsolation()
    
    #=====prepare biophysical parameters at overpass time======================
    sceneDir = os.path.join(landsatDataBase,'albedo',scene)
    outFN = os.path.join(sceneDir,'%s_albedo.tiff' % sceneID) 
    if not os.path.exists(outFN):
        if not os.path.exists(sceneDir):
            os.makedirs(sceneDir)
        print 'processing : albedo...' 
        a = Landsat(fn,LCpath)
        a.getAlbedo()

    
    sceneDir = os.path.join(landsatDataBase,'LC',scene)
    outFN = os.path.join(sceneDir,'%s_LC.tiff' % sceneID)
    if not os.path.exists(outFN):
        if not os.path.exists(sceneDir):
            os.makedirs(sceneDir)
        a = Landsat(fn,LCpath)
        a.getLC(landcover)

def main():    
    # Get time and location from user
    parser = argparse.ArgumentParser()
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("lon", type=float, help="longitude")
    parser.add_argument("isUSA", type=float, help="USA=1, non-USA=0")
    parser.add_argument("start_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("end_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("ET_dir", type=str, help="ALEXI ET directory")
    parser.add_argument("LC_dir", type=str, help="Landcover directory")
    parser.add_argument("cloud", type=int, help="cloud cover")
    parser.add_argument("collection", type=int,nargs='?', default=1)
    parser.add_argument('-s','--sat', nargs='?',type=int, default=8, help='which landsat to search or download, i.e. Landsat 8 = 8')
    args = parser.parse_args()
    loc = [args.lat,args.lon] 
    isUSA = args.isUSA
    start_date = args.start_date
    end_date = args.end_date
    ET_dir = args.ET_dir
    LC_dir = args.LC_dir
    cloud = args.cloud
    sat = args.sat
    collection = args.collection
      
     # =====earthData credentials==============================================
    earth_user = str(getpass.getpass(prompt="earth login username:"))
    if keyring.get_password("nasa",earth_user)==None:
        earth_pass = str(getpass.getpass(prompt="earth login password:"))
        keyring.set_password("nasa",earth_user,earth_pass)
    else:
        earth_pass = str(keyring.get_password("nasa",earth_user)) 
        
    # =====USGS credentials====================================================
     # need to get this from pop up
    usgs_user = str(getpass.getpass(prompt="usgs username:"))
    if keyring.get_password("usgs",usgs_user)==None:
        usgs_pass = str(getpass.getpass(prompt="usgs password:"))
        keyring.set_password("usgs",usgs_user,usgs_pass)
    else:
        usgs_pass = str(keyring.get_password("usgs",usgs_user)) 
        
    session = (earth_user, earth_pass)
#    collection = 1

    #===process Landsat LAI====================================================
    print("processing LAI...")
#    processlai.get_LAI(loc,start_date,end_date,usgs_user,
#                       usgs_pass,earth_user,earth_pass,cloud)
    
    processlai.get_LAI(loc,start_date,end_date,earth_user,
                       earth_pass,cloud,sat,cacheDir)
    
    #===process met,alexi and misc landsat data================================
    print("processing MET,ALEXI and misc landsat data ...")
#    landsatTemp = os.path.join(landsatSR,'temp')
#    fileList = glob.glob(os.path.join(landsatTemp,"*_MTL.txt"))
    available = 'Y'
    Downloaded_df = getlandsatdata.search(loc[0],loc[1],start_date,end_date,cloud,available,cacheDir,sat)
    fileList = Downloaded_df.local_file_path
    for fn in fileList:
        prepare_data(fn,session,isUSA,LC_dir,ET_dir)
    
    #===process Landsat LST====================================================
    print("processing LST...")
#    processlst.get_lst(earth_user,earth_pass)
    processlst.get_lst(loc,start_date,end_date,earth_user,earth_pass,cloud,sat,cacheDir)
    
if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)  
    
    
    
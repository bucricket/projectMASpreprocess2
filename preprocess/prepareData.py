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
import pycurl
from .utils import folders
from .processData import ALEXI,MET,Landsat
from .landsatTools import landsat_metadata
from processlai import processlai
from processlst import processlst
from getlandsatdata import getlandsatdata
import pandas as pd
import sqlite3
from os import walk
import shutil


#===set golbal paths===========================================================
base = os.getcwd()
cacheDir = os.path.abspath(os.path.join(base,os.pardir,"SATELLITE_DATA"))
model_cache = os.path.abspath(os.path.join(base,os.pardir,"MODEL_DATA"))
Folders = folders(base)  
ALEXIbase = Folders['ALEXIbase']
metBase = Folders['metBase']
landsatDataBase = Folders['landsatDataBase']

def moveFiles(top_path,dst_path,ext):           
    for root, dirs, files in os.walk(top_path):
        for file in files:
            if file.endswith(ext):            
                 srccpy = os.path.join(root, file)
                 dstcpy = os.path.join(dst_path,file)
                 if not os.path.exists(dstcpy):
                     shutil.copy(srccpy,dstcpy)

def prepare_data(fn,session,isUSA,LCpath):
    

    LCpath = LCpath
#    ETpath = ETpath
    meta = landsat_metadata(fn)
    if meta.SPACECRAFT_ID == 'LANDSAT_7':        
        sat = 7
    else:
        sat = 8
    sceneID = meta.LANDSAT_SCENE_ID
    productID = meta.LANDSAT_PRODUCT_ID
    landsatCacheDir = os.path.join(cacheDir,"LANDSAT")
    satscene_path = os.sep.join(fn.split(os.sep)[:-2])
#    scene = sceneID[3:9]
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
    output_df = getlandsatdata.searchProduct(productID,landsatCacheDir,sat)
        
    #=====pick Landcover map===================================================
    if isUSA ==1:
        landcover = 'NLCD'
    else:
        landcover = 'GlobeLand30'

    
    #=====prepare the ETd data=================================================    
#    sceneDir = os.path.join(ALEXIbase,'%s' % scene) 
    sceneDir = os.path.join(satscene_path,'ET','400m')   
    if not os.path.exists(sceneDir):
        os.mkdir(sceneDir)     
    
#    outFN = os.path.join(sceneDir,'%s_alexiETSub.tiff' % sceneID) 
    outFN = os.path.join(sceneDir,'%s_alexiET.tiff' % sceneID) 
    if not os.path.exists(outFN):
        print 'get->ALEXI ET...'
        a = ALEXI(fn)
        a.getALEXIdata(ALEXIgeodict,isUSA)
        processlai.updateLandsatProductsDB(output_df,outFN,landsatCacheDir,'ALEXI_ET')
    
    #=====prepare MET data=====================================================    
#    sceneDir = os.path.join(metBase,'%s' % scene)
    sceneDir = os.path.join(satscene_path,'MET')
    if not os.path.exists(sceneDir):
        os.mkdir(sceneDir) 
    
#    outFN = os.path.join(sceneDir,'%s_pSub.tiff' % sceneID) 
    outFN = os.path.join(sceneDir,'%s_p.tiff' % sceneID) 
    outFNu = os.path.join(sceneDir,'%s_u.tiff' % sceneID)
    outFNta = os.path.join(sceneDir,'%s_Ta.tiff' % sceneID)
    outFNq2 = os.path.join(sceneDir,'%s_q2.tiff' % sceneID)
    if not os.path.exists(outFN):
        print 'get->MET data...'
        a = MET(fn,session)
        a.getCFSR()
        processlai.updateLandsatProductsDB(output_df,outFN,landsatCacheDir,'SFC_PRESS')
        processlai.updateLandsatProductsDB(output_df,outFNu,landsatCacheDir,'WIND')
        processlai.updateLandsatProductsDB(output_df,outFNta,landsatCacheDir,'TA')
        processlai.updateLandsatProductsDB(output_df,outFNq2,landsatCacheDir,'Q2')
    
    #====prepare insolation====================================================
    sceneDir = os.path.join(satscene_path,'INSOL')
    if not os.path.exists(sceneDir):
        os.mkdir(sceneDir) 
#    outFN = os.path.join(sceneDir,'%s_Insol1Sub.tiff' % sceneID)
    outFN = os.path.join(sceneDir,'%s_Insol1.tiff' % sceneID)
    outFN24 = os.path.join(sceneDir,'%s_Insol24.tiff' % sceneID)
    if not os.path.exists(outFN):
        a = MET(fn,session)
#        a.getInsolation()
        a.getGSIP()
        processlai.updateLandsatProductsDB(output_df,outFN,landsatCacheDir,'INSOL1')
        processlai.updateLandsatProductsDB(output_df,outFN24,landsatCacheDir,'INSOL24')
    
    #=====prepare biophysical parameters at overpass time======================
#    sceneDir = os.path.join(landsatDataBase,'albedo',scene)
    sceneDir = os.path.join(satscene_path,'ALBEDO')
    if not os.path.exists(sceneDir):
        os.mkdir(sceneDir) 
    outFN = os.path.join(sceneDir,'%s_albedo.tiff' % sceneID) 
    if not os.path.exists(outFN):
        print 'processing : albedo...' 
        a = Landsat(fn,LCpath)
        a.getAlbedo()
        processlai.updateLandsatProductsDB(output_df,outFN,landsatCacheDir,'ALBEDO')

    
#    sceneDir = os.path.join(landsatDataBase,'LC',scene)
    sceneDir = os.path.join(satscene_path,'LC')
    if not os.path.exists(sceneDir):
        os.mkdir(sceneDir) 
    outFN = os.path.join(sceneDir,'%s_LC.tiff' % sceneID)
    if not os.path.exists(outFN):
        a = Landsat(fn,LCpath)
        a.getLC(landcover)
        processlai.updateLandsatProductsDB(output_df,outFN,landsatCacheDir,'LC')

def main():    
    # Get time and location from user
    parser = argparse.ArgumentParser()
    parser.add_argument("lat", type=float, help="latitude")
    parser.add_argument("lon", type=float, help="longitude")
    parser.add_argument("isUSA", type=float, help="USA=1, non-USA=0")
    parser.add_argument("start_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("end_date", type=str, help="Start date yyyy-mm-dd")
    parser.add_argument("LC_dir", type=str, help="Landcover directory")
    parser.add_argument("cloud", type=int, help="cloud cover")
    dst_path = os.path.join(model_cache,"ALEXI")
    parser.add_argument("--ET_dir", type=str,default=dst_path,help="ALEXI ET top directory")
    dst_path = os.path.join(cacheDir,"GSIP")
    parser.add_argument("--Insol_dir", type=str,default=dst_path,help="Insolation top directory")
    parser.add_argument("collection", type=int,nargs='?', default=1)
    parser.add_argument('-s','--sat', nargs='?',type=int, default=8, help='which landsat to search or download, i.e. Landsat 8 = 8')
    args = parser.parse_args()
    loc = [args.lat,args.lon] 
    isUSA = args.isUSA
    start_date = args.start_date
    end_date = args.end_date
    ET_dir = args.ET_dir
    Insol_dir = args.Insol_dir
    LC_dir = args.LC_dir
    cloud = args.cloud
    sat = args.sat
    collection = args.collection
    #=======copy ALEXI ET and Insol files to central cache=====================
    ext = ".tif"
    dst_path = os.path.join(model_cache,"ALEXI")
    if not os.path.exists(dst_path):
            os.mkdir(dst_path)
    moveFiles(ET_dir,dst_path,ext)
    ext = ".gz"
    dst_path = os.path.join(cacheDir,"GSIP")
    if not os.path.exists(dst_path):
            os.mkdir(dst_path)
    moveFiles(Insol_dir,dst_path,ext)
     
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
    processlai.get_LAI(loc,start_date,end_date,earth_user,
                       earth_pass,cloud,sat,cacheDir)
    
    #===process met,alexi and misc landsat data================================
    print("processing MET,ALEXI and misc landsat data ...")
    landsatCacheDir = os.path.join(cacheDir,"LANDSAT")
    db_fn = os.path.join(landsatCacheDir,"landsat_products.db")
    available = 'Y'
    product = 'ALBEDO'
    search_df = getlandsatdata.search(loc[0],loc[1],start_date,end_date,cloud,available,landsatCacheDir,sat)
    productIDs = search_df.LANDSAT_PRODUCT_ID
#    fileList = search_df.local_file_path
    #====check what products are processed against what Landsat data is available===
    if os.path.exists(db_fn):
        conn = sqlite3.connect( db_fn )
        res = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = res.fetchall()[0]
        if (product in tables):
            processedProductIDs = processlst.searchLandsatProductsDB(loc[0],loc[1],start_date,end_date,product,landsatCacheDir)
            df1 = processedProductIDs[["LANDSAT_PRODUCT_ID"]]
            merged = df1.merge(pd.DataFrame(productIDs), indicator=True, how='outer')
            df3 = merged[merged['_merge'] != 'both' ]
            productIDs = df3[["LANDSAT_PRODUCT_ID"]].LANDSAT_PRODUCT_ID
            if len(productIDs)>0:
                output_df = pd.DataFrame()
                for productID in productIDs:
                    output_df = output_df.append(getlandsatdata.searchProduct(productID,cacheDir,sat),ignore_index=True)
                productIDs = output_df.LANDSAT_PRODUCT_ID

    for productID in productIDs:
        out_df = getlandsatdata.searchProduct(productID,landsatCacheDir,sat)
        fn = os.path.join(out_df.local_file_path[0],productID+"_MTL.txt")
        print fn
        prepare_data(fn,session,isUSA,LC_dir)
    
    #===process Landsat LST====================================================
    print("processing LST...")
#    processlst.get_lst(earth_user,earth_pass)
    processlst.get_lst(loc,start_date,end_date,earth_user,earth_pass,cloud,sat,cacheDir)
    
if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)  
    
    
    
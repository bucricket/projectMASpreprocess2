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
import numpy as np
import wget
from .utils import folders
from .processData import ALEXI, MET, Landsat
from .landsatTools import landsat_metadata
from processlai import processlai
from processlst import processlst
from getlandsatdata import getlandsatdata
import pandas as pd
import sqlite3
from os import walk
import shutil
import datetime
import glob

# ===set golbal paths===========================================================
base = os.getcwd()
cacheDir = os.path.abspath(os.path.join(base, "SATELLITE_DATA"))
model_cache = os.path.abspath(os.path.join(base, "MODEL_DATA"))
Folders = folders(base)
ALEXIbase = Folders['ALEXIbase']
metBase = Folders['metBase']
landsatDataBase = Folders['landsatDataBase']


def moveFiles(top_path, dst_path, ext):
    for root, dirs, files in os.walk(top_path):
        for file in files:
            if file.endswith(ext):
                srccpy = os.path.join(root, file)
                dstcpy = os.path.join(dst_path, file)
                if not os.path.exists(dstcpy):
                    shutil.copy(srccpy, dstcpy)


def searchLandsatSceneID(sceneID, db_path, sat):
    columns = ['acquisitionDate', 'acquisitionDate', 'upperLeftCornerLatitude', 'upperLeftCornerLongitude',
               'lowerRightCornerLatitude', 'lowerRightCornerLongitude', 'cloudCover']
    if sat == 7:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_ETM_C1.csv'
        db_name = os.path.join(db_path, 'LANDSAT_ETM_C1.db')
    else:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_8_C1.csv'
        db_name = os.path.join(db_path, 'LANDSAT_8_C1.db')

    fn = os.path.join(db_path, metadataUrl.split(os.sep)[-1])
    if not os.path.exists(db_name):
        if not os.path.exists(fn):
            wget.download(metadataUrl, out=fn)
        conn = sqlite3.connect(db_name)
        orig_df = pd.read_csv(fn, usecols=columns)
        orig_df['sr'] = pd.Series(np.tile('N', len(orig_df)))
        orig_df['bt'] = pd.Series(np.tile('N', len(orig_df)))
        orig_df['local_file_path'] = ''
        orig_df.to_sql("raw_data", conn, if_exists="replace", index=False)
        conn.close()

    conn = sqlite3.connect(db_name)
    output = pd.read_sql_query("SELECT * from raw_data WHERE (sceneID == '%s')" % sceneID, conn)
    conn.close()
    return output


def updateLandsatProductsDB(cacheDir, product, fn):
    #
    #    if product == 'LST':
    #        product_name = 'lstSharp'
    #    elif product == 'LAI':
    #        product_name = 'lai'
    #    elif product == 'ALBEDO':
    #        product_name = 'albedo'
    #    elif product == 'CF_MASK':
    #        product_name = 'Mask'
    #    elif product == 'INSOL1':
    #        product_name = 'Insol1'
    #    elif product == 'INSOL24':
    #        product_name = 'Insol24'
    #    elif product == 'LC':
    #        product_name = 'LC'
    #    elif product == 'SFC_PRESS':
    #        product_name = 'p'
    #    elif product == 'WIND':
    #        product_name = 'u'
    #    elif product == 'TA':
    #        product_name = 'Ta'
    #    elif product == 'Q2':
    #        product_name = 'q2'
    #    elif product == 'NDVI':
    #        product_name = 'ndvi'
    #    elif product == 'ALEXI_ET':
    #        product_name = 'alexiET'

    #    filenames =[]
    #    paths =[]
    #    for dirpath, dirnames, fns in os.walk(topDir):
    #        try:
    #            for filename in [f for f in fns if ((f.split("_")[1].split(".")[0] == "%s" % product_name) and f.endswith(".tiff"))]:
    #                filenames.append(filename)
    #                paths.append(os.path.join(dirpath,filename))
    #        except:
    #            pass
    #
    #    landsatDB = pd.DataFrame()
    #    for fn in filenames:
    sceneID = fn.split("_")[0]
    sat = sceneID[2]
    landsatDB = searchLandsatSceneID(sceneID, cacheDir, sat)

    if not len(fn) == 0:
        db_fn = os.path.join(cacheDir, "landsat_products.db")

        date = landsatDB.acquisitionDate
        ullat = landsatDB.upperLeftCornerLatitude
        ullon = landsatDB.upperLeftCornerLongitude
        lllat = landsatDB.lowerRightCornerLatitude
        lllon = landsatDB.lowerRightCornerLongitude
        productIDs = landsatDB.LANDSAT_PRODUCT_ID

        if not os.path.exists(db_fn):
            conn = sqlite3.connect(db_fn)
            landsat_dict = {"acquisitionDate": date, "upperLeftCornerLatitude": ullat,
                            "upperLeftCornerLongitude": ullon,
                            "lowerRightCornerLatitude": lllat,
                            "lowerRightCornerLongitude": lllon,
                            "LANDSAT_PRODUCT_ID": productIDs, "filename": fn}
            landsat_df = pd.DataFrame.from_dict(landsat_dict)
            landsat_df.to_sql("%s" % product, conn, if_exists="replace", index=False)
            conn.close()
        else:
            conn = sqlite3.connect(db_fn)
            res = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = res.fetchall()[0]
            if (product in tables):
                orig_df = pd.read_sql_query("SELECT * from %s" % product, conn)
            else:
                orig_df = pd.DataFrame()

            landsat_dict = {"acquisitionDate": date, "upperLeftCornerLatitude": ullat,
                            "upperLeftCornerLongitude": ullon,
                            "lowerRightCornerLatitude": lllat,
                            "lowerRightCornerLongitude": lllon,
                            "LANDSAT_PRODUCT_ID": productIDs, "filename": fn}
            landsat_df = pd.DataFrame.from_dict(landsat_dict)
            orig_df = orig_df.append(landsat_df, ignore_index=True)
            orig_df = orig_df.drop_duplicates(keep='last')
            orig_df.to_sql("%s" % product, conn, if_exists="replace", index=False)
            conn.close()


def search(lat, lon, start_date, end_date, cloud, cacheDir, sat):
    columns = ['acquisitionDate', 'acquisitionDate', 'upperLeftCornerLatitude', 'upperLeftCornerLongitude',
               'lowerRightCornerLatitude', 'lowerRightCornerLongitude', 'cloudCover', 'sensor', 'LANDSAT_PRODUCT_ID']
    end = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    # this is a landsat-util work around when it fails
    if sat == 7:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_ETM_C1.csv'
    else:
        metadataUrl = 'https://landsat.usgs.gov/landsat/metadata_service/bulk_metadata_files/LANDSAT_8_C1.csv'

    fn = os.path.join(cacheDir, metadataUrl.split(os.sep)[-1])
    # looking to see if metadata CSV is available and if its up to the date needed
    if os.path.exists(fn):
        d = datetime.datetime.fromtimestamp(os.path.getmtime(fn))
        if (end.year > d.year) and (end.month > d.month) and (end.day > d.day):
            wget.download(metadataUrl, out=fn)
            df = pd.read_csv(fn, usecols=columns)
            df.to_csv(fn)
        df = pd.read_csv(fn)
        index = ((df.acquisitionDate >= start_date) & (df.acquisitionDate < end_date) & (
                df.upperLeftCornerLatitude > lat) & (df.upperLeftCornerLongitude < lon) & (
                         df.lowerRightCornerLatitude < lat) & (df.lowerRightCornerLongitude > lon) & (
                         df.cloudCover <= cloud) & (df.sensor == 'OLI_TIRS'))
        df = df[index]

    else:
        wget.download(metadataUrl, out=fn)
        df = pd.read_csv(fn, usecols=columns)
        df.to_csv(fn)
        index = ((df.acquisitionDate >= start_date) & (df.acquisitionDate < end_date) & (
                df.upperLeftCornerLatitude > lat) & (df.upperLeftCornerLongitude < lon) & (
                         df.lowerRightCornerLatitude < lat) & (df.lowerRightCornerLongitude > lon) & (
                         df.cloudCover <= cloud) & (df.sensor == 'OLI_TIRS'))
        df = df[index]

    return df


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def find_already_downloaded(df, cache_dir):
    usgs_available = list(df.LANDSAT_PRODUCT_ID.values)
    # find sat
    sat = usgs_available[0].split("_")[0][-1]
    # find scenes
    scenes = [x.split("_")[2] for x in usgs_available]
    scenes = list(set(scenes))
    available_list = []
    for scene in scenes:
        path_to_search = os.path.join(cache_dir, 'L%s/%s/RAW_DATA/*MTL*' % (sat, scene))
        available = [os.path.basename(x) for x in
                     glob.glob(path_to_search)]
        available = [x[:-8] for x in available]
        available_list = available_list + available
    return intersection(usgs_available, available_list)


def find_not_processed(downloaded, cache_dir):
    """finds the files that are downloaded but still need to process ALBEDO data and thus the rest of the inputs"""
    # find sat
    sat = downloaded[0].split("_")[0][-1]
    # find scenes
    scenes = [x.split("_")[2] for x in downloaded]
    scenes = list(set(scenes))
    available_list = []
    for scene in scenes:
        path_to_search = os.path.join(cache_dir, 'L%s/%s/ALBEDO/*_l.tif' % (sat, scene))
        available = [os.path.basename(x) for x in
                     glob.glob(path_to_search)]
        available = [x[:-8] for x in available]
        available_list = available_list + available
    for x in available_list:
        if x in downloaded:
            downloaded.remove(x)
    return downloaded


def prepare_data(fn, session, isUSA, LCpath, insolDataset):
    LCpath = LCpath
    meta = landsat_metadata(fn)
    if meta.SPACECRAFT_ID == 'LANDSAT_7':
        sat = 7
    else:
        sat = 8
    sceneID = meta.LANDSAT_SCENE_ID
    satscene_path = os.sep.join(fn.split(os.sep)[:-2])
    isUSA = isUSA
    ALEXI_ulLon = 0.0
    ALEXI_ulLat = 0.0
    # --------------------
    ALEXILatRes = 0.004
    ALEXILonRes = 0.004
    ALEXIshape = [3750, 3750]
    ALEXIgeodict = {'ALEXI_ulLat': ALEXI_ulLat, 'ALEXI_ulLon': ALEXI_ulLon,
                    'ALEXI_LatRes': ALEXILatRes, 'ALEXI_LonRes': ALEXILonRes,
                    'ALEXIshape': ALEXIshape}
    # =====pick Landcover map===================================================
    if isUSA == 1:
        landcover = 'NLCD'
    else:
        landcover = 'GlobeLand30'

    # =====prepare the ETd data=================================================
    #    sceneDir = os.path.join(ALEXIbase,'%s' % scene)
    sceneDir = os.path.join(satscene_path, 'ET', '400m')
    if not os.path.exists(sceneDir):
        os.makedirs(sceneDir)
    outFN = os.path.join(sceneDir, '%s_alexiET.tiff' % sceneID)
    if not os.path.exists(outFN):
        print 'get->ALEXI ET...'
        a = ALEXI(fn)
        a.getALEXIdata(ALEXIgeodict, isUSA)
    # =====prepare MET data=====================================================
    sceneDir = os.path.join(satscene_path, 'MET')
    if not os.path.exists(sceneDir):
        os.makedirs(sceneDir)
    outFN = os.path.join(sceneDir, '%s_p.tiff' % sceneID)
    if not os.path.exists(outFN):
        print 'get->MET data...'
        a = MET(fn, session)
        a.getCFSR()

    # ====prepare insolation====================================================
    sceneDir = os.path.join(satscene_path, 'INSOL')
    if not os.path.exists(sceneDir):
        os.makedirs(sceneDir)
    outFN = os.path.join(sceneDir, '%s_Insol1.tiff' % sceneID)
    if not os.path.exists(outFN):
        a = MET(fn, session)
        #        a.getInsolation()
        if insolDataset == 'GSIP':
            a.getGSIP()
        elif insolDataset == 'CERES':
            a.getCERESinsol()
        else:
            print("Using MERRA for Insolation")
            a.getInsolation()

    # =====prepare biophysical parameters at overpass time======================
    sceneDir = os.path.join(satscene_path, 'ALBEDO')
    if not os.path.exists(sceneDir):
        os.makedirs(sceneDir)
    outFN = os.path.join(sceneDir, '%s_albedo.tiff' % sceneID)
    if not os.path.exists(outFN):
        print 'processing : albedo...'
        a = Landsat(fn, LCpath)
        a.getAlbedo()

    sceneDir = os.path.join(satscene_path, 'LC')
    if not os.path.exists(sceneDir):
        os.makedirs(sceneDir)
    outFN = os.path.join(sceneDir, '%s_LC.tiff' % sceneID)
    if not os.path.exists(outFN):
        a = Landsat(fn, LCpath)
        a.getLC(landcover)


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
    parser.add_argument("insolDataset", type=str, help="insolation dataset: CERES or GSIP")
    dst_path = os.path.join(model_cache, "ALEXI")
    parser.add_argument("--ET_dir", type=str, default=dst_path, help="ALEXI ET top directory")
    dst_path = os.path.join(cacheDir, "MERRA")
    parser.add_argument("--Insol_dir", type=str, default=dst_path, help="Insolation top directory")
    parser.add_argument('-s', '--sat', nargs='?', type=int, default=8,
                        help='which landsat to search or download, i.e. Landsat 8 = 8')
    args = parser.parse_args()
    loc = [args.lat, args.lon]
    isUSA = args.isUSA
    start_date = args.start_date
    end_date = args.end_date
    ET_dir = args.ET_dir
    Insol_dir = args.Insol_dir
    LC_dir = args.LC_dir
    cloud = args.cloud
    sat = args.sat
    insolDataset = args.insolDataset

    # =======copy ALEXI ET and Insol files to central cache=====================
    ext = ".tif"
    dst_path = os.path.join(model_cache, "ALEXI")
    if not os.path.exists(dst_path):
        os.mkdir(dst_path)
    moveFiles(ET_dir, dst_path, ext)
    if insolDataset == 'GSIP':
        ext = ".gz"
        dst_path = os.path.join(cacheDir, "GSIP")
        if not os.path.exists(dst_path):
            os.makedirs(dst_path)
        moveFiles(Insol_dir, dst_path, ext)
    elif insolDataset == 'CERES':
        ext = ".nc"
        dst_path = os.path.join(cacheDir, "CERES")
        if not os.path.exists(dst_path):
            os.makedirs(dst_path)
        moveFiles(Insol_dir, dst_path, ext)
    else:
        ext = ".tif"

    # =====earthData credentials==============================================
    earth_user = str(getpass.getpass(prompt="earth login username:"))
    if keyring.get_password("nasa", earth_user) == None:
        earth_pass = str(getpass.getpass(prompt="earth login password:"))
        keyring.set_password("nasa", earth_user, earth_pass)
    else:
        earth_pass = str(keyring.get_password("nasa", earth_user))

        # =====USGS credentials====================================================
    # need to get this from pop up
    usgs_user = str(getpass.getpass(prompt="usgs username:"))
    if keyring.get_password("usgs", usgs_user) == None:
        usgs_pass = str(getpass.getpass(prompt="usgs password:"))
        keyring.set_password("usgs", usgs_user, usgs_pass)
    else:
        usgs_pass = str(keyring.get_password("usgs", usgs_user))

    session = (earth_user, earth_pass)
    #    collection = 1

    # ===process Landsat LAI====================================================
    print("processing LAI...")
    processlai.get_LAI(loc, start_date, end_date, earth_user,
                       earth_pass, cloud, sat, cacheDir)

    # ===process met,alexi and misc landsat data================================
    print("processing MET,ALEXI and misc landsat data ...")
    landsatCacheDir = os.path.join(cacheDir, "LANDSAT")
    output_df = search(loc[0], loc[1], start_date, end_date, cloud, landsatCacheDir, sat)
    downloaded = find_already_downloaded(output_df, landsatCacheDir)
    productIDs = find_not_processed(downloaded, landsatCacheDir)

    for productID in productIDs:
        sat_str = productID.split("_")[0][-1]
        scene = productID.split("_")[2]
        path = os.path.join(landsatCacheDir, 'L%s/%s/RAW_DATA/' % (sat_str, scene))
        fn = os.path.join(path, productID + "_MTL.txt")
        print fn
        prepare_data(fn, session, isUSA, LC_dir, insolDataset)

    # ===process Landsat LST====================================================
    print("processing LST...")
    #    processlst.get_lst(earth_user,earth_pass)
    processlst.get_lst(loc, start_date, end_date, earth_user, earth_pass, cloud, sat, cacheDir)


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, pycurl.error):
        exit('Received Ctrl + C... Exiting! Bye.', 1)

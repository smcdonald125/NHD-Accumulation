"""
python /media/imagery/usgs_sc/smcdonald/NHDv2_Agg/ForGitHub/version2/accumulateNHDPlus.py
Purpose:Calculate raster statistics at the NHD+ catchment level and accumulate the data. The input raster can be continuous data,
        where the following statistics can be calculated: MAX, MIN, MEAN, MEDIAN, SUM. The raster can also be classed data, where the
        total area per class is calculated. The catchment statistics and accumulated data can be written as a table to a CSV and
        can be written as a shapefile, joined to the NHD+ catchments. The upstream networks are built using a routing database, produced
        by Mike Wieczorek. The routing database can be reduced to the VPU zone level, and can be read in as its original form (SASS database)
        or a CSV file. 
Data:   NHD+ catchments as shapefile - Shapefile to act as zones - must have FEATUREID and AreaSqKM fields
        Routing Database as SASS database or CSV- Must have FL_ComID, cfromnode, ctonode, nhdplusreg fields
        Data to accumulate:
            Raster File (tif) - user must specify if the data is continuous or classed 
              OR
            Tabular data (shp or csv) - tabular data containing values for each catchment
Note: The user is responsible for ensuring that the catchments and the raster are in the same projection.
Authors:  Sarah McDonald, Independent Contractor, USGS
Contact: smcdonald@chesapeakebay.net
"""

#Import Libs
import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.mask
import fiona
from osgeo import gdal
import os
from sas7bdat import SAS7BDAT
import time

import snapped_raster_stats as srs

####################################################################################################
####################       READ INPUT DATA       ###################################################
####################################################################################################
def readCatchmentGeoms(NHDpath, catchID):
    """
    Method: readCatchmentGeoms()
    Purpose: Read NHD+ catchment geometries into a dictionary
    Params: NHDpath - path to NHD+ catchments (shapefile)
            catchID - the name of the column in the shapefile to use as the unique ID
    Returns: shapes - dictionary of fiona geometries for NHD+ catchments
    """
    print("Reading catchment geometries...")
    if os.path.isfile(NHDpath):
        gdf = gpd.read_file(NHDpath)
        print(list(gdf))
        crs = gdf.crs
        del gdf
        with fiona.open(NHDpath, "r") as geoms:
            shapes = {feature['properties'][catchID]:feature['geometry'] for feature in geoms}
        return shapes, crs
    else:
        print("Invalid File Path: ", NHDpath)
        print("Could not read catchments shapefile")
        return {}

def readRouting(sassDB, VPUzone):
    """
    Method: readRouting()
    Purpose: Read routing database into a dataframe and return it.
    Params: sassDB - path to routing database (SASS database or CSV)
            VPUzone - VPU zone (from NHD) that the catchments exist in. Used to shrink 
                    routing database from CONUS scale.
    Returns: db - pandas dataframe of routing database
    """
    print("Reading routing database...")
    if os.path.isfile(sassDB):
        if sassDB[-8:] == "sas7bdat":
            with SAS7BDAT (sassDB, skip_header=False) as r:
                # cols = r.header.Name #list of all columns
                db = r.to_data_frame()
                print(list(db))
                db = db[['FL_ComID', 'cfromnode', 'ctonode', 'nhdplusreg']]
                db = db[db.nhdplusreg == VPUzone]
                db.to_csv("/data/gis2/sm/NHD_Test/Version3_MWDB/Output/SassDB_02.csv", index=False)
                db = db.fillna(-100)
                db["cfromnode"] = db.cfromnode.astype(int)
                db["ctonode"] = db.ctonode.astype(int)
        elif sassDB[-3:] == "csv":
            db = pd.read_csv(sassDB)
            db = db.dropna(axis='index')
            db['FL_ComID'] = db['FL_ComID'].apply(int)
        else:
            print("Invalid routing database path\n\tMust be .sas7bdat file or .csv file")
            return pd.DataFrame()
        return db
    else:
        print("Invalid File Path: ", sassDB)
        return pd.DataFrame()

def readRas(rasPath):
    """
    Method: readRas()
    Purpose: Verify path and get cell size and no data for the raster.
    Params:  rasPath - path to raster
    Returns: cellSize - cell size of the raster
             noData - noData value of the raster
    """
    if os.path.isfile(rasPath):
        with rio.open(rasPath, 'r') as src:
            noData = src.nodatavals
            noData = noData[0]
            t = src.transform
            cellSize = t.a

        return cellSize, noData     
    else:
        print("Invalid file path: ", rasPath)
        sys.exit()

def readTabular(tabPath, idName, cols, PLUS):
    """
    Method: readTabular()
    Purpose: create data frame with catchment statistics from an existing table. User specifies column names
            to accumulate and the column name denoting the unique catchment IDs.
    Params: tabPath - path to the raster data
            idName - string of column name containing the unique COMIDs
            cols - list of columns names to accumulate
    Returns: sumDict - dictionary where the catchment ID is the key and the data is a list of catchment statistics and
                    accumulated statistics
    """
    print("Reading statistics from table...")
    if os.path.isfile(tabPath):
        if tabPath[-3:] == 'csv':
            df = pd.read_csv(tabPath)
        elif tabPath[-3:] == 'shp':
            df = gpd.read_file(tabPath)
        else:
            print("Invalid file type: ", tabPath)
            print("Tabular data must be CSV or shp")
            print(tabPath[-3:])
            return {}
        #remove unecessary columns and convert to dictionary
        df = df[[idName] + cols]
        df[idName] = df[idName].astype(int)
        ids = list(set(list(PLUS['FL_ComID'])))
        orig_len = len(df)
        df = df[df[idName].isin(ids)]
        if len(df) < orig_len:
            dif = orig_len - len(df)
            print(dif, "Records removed - outside watershed")
        df = df.set_index(idName)
        df = df[cols]
        df = df.T #transpose data so that each catchment ID is a column
        c = df.columns.values
        df = df.set_index(c[0])
        df[c[0]] = df.index
        sumDict = df.to_dict('list')
        return sumDict
    else:
        print("Cannot read statistics table")
        print("Invalid file path: ", tabPath)
        return {}

####################################################################################################
####################       ACCUMULATE DATA       ###################################################
####################################################################################################
def accumulate_mp(rasType, PLUS, dirLinks, network, sumDict, convFact, direction):
    """
    Method: accumulate_mp()
    Purpose: Divide the catchments into chunks and multiprocess the accumulation of raster statistics
    Params: rasType - 0 or 1 specifying if raster is continuous or classed
            PLUS - dataframe of the routing database
            dirLinks - dictionary of direct network links
            network - dictionary storing full network in specified direction
            sumDict - dictionary storing catchment stats
                    Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                    Classed data: dataframe of class value and pixel count
            convFact - the conversion factor specified by the user
            direction - string upstream or downstream
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
            network - dictionary storing full network in specified direction
    """
    cpus_minus_1 = mp.cpu_count() - 2
    if cpus_minus_1 == 0:
        cpus_minus_1 = 1

    allCatch = list(sumDict)
    
    chunk_iterator = []
    for i in range(len(allCatch)):
        gdf_args = allCatch[i], PLUS, dirLinks, network, sumDict, convFact, direction, rasType
        chunk_iterator.append(gdf_args)

    pool = mp.Pool(processes=cpus_minus_1)

    results = dict()
    totalResults = 0
    acc_res = pool.map(accumulate, chunk_iterator)
    for result, upresult in acc_res: #list of dictionaries
        totalResults += len(result)
        results.update(result)
        if totalResults != len(results):
            print("Dict was not updated")
        network.update(upresult)
    pool.close()

    return results, network

def accumulate(args): 
    """
    Method: accumulate()
    Purpose: accumulate the raster statistics
    Params: allCatch - list of catchment IDs to loop through and accumulate
            PLUS - dataframe of the routing database
            dirLinks - dictionary of direct network links
            network - dictionary storing full network in specified direction
            sumDict - dictionary storing catchment stats
                    Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                    Classed data: dataframe of class value and pixel count
            convFact - the conversion factor specified by the user
            direction - string upstream or downstream
            rasType - 0 or 1 specifying if raster is continuous or classed
            toFromNames - list with column names; index 0 is downstream connection name, index 1 is current catchment name
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
    """
    allCatch, PLUS, dirLinks, network, sumDict, convFact, direction, rasType = args
    #Accumulate statistic for each catchment
    allSum = list(sumDict)
    newSum, newUp = {}, {}
    for s in allCatch:
        if s not in network:
            network[s] = generateNetwork(dirLinks, s, PLUS, direction) #all network connections
            newUp[s] = network[s].copy()
        upCatch = list((set(network[s])&set(allSum))) #all network catchments that are in sumDict
        newSum[s] = accumulate_values(sumDict, upCatch, convFact, s, rasType) #list of catchment values (original scale) and accumulated values (converted scale)
    return newSum, newUp #return all data in dict

def accumulate_values(sumDict, toSum, convFact, curID, rasType): 
    """
    Method: accumulate_values()
    Purpose: Sum each statistic value for all catchments in the network.
    Params: sumDict - dictionary containing catchment statistics
            toSum - list of network catchment IDs
            convFact - the conversion factor to be applied to the accumulated values 
            curID - the "current" catchment ID
            rasType - 0 or 1 specifying if raster is continuous or classed
    Returns: allStats - list of catchment statistics and accumulated statistics for the current catchment
                        OR
             df - dataframe of accumulated catchment statistics
    """
    if rasType == 0: # continuous raster
        sums = [0] * 2
        for cID in toSum:
            if cID in sumDict: #check that network catchment is in the dict
                tmp = sumDict[cID] #get list of all stats for upstream catchment
                for idx in range(2): #for stat and area - sum
                    sums[idx] += tmp[idx]
        for i in range(len(sums)):
            sums[i] = sums[i] * convFact #convert data
        allStats = sumDict[curID] + sums #create full list in the order "catchment statistics for each class and accumulated stats for each class"
        return allStats
    else: # classed data
        tmp = []
        for cID in toSum:
            if cID in sumDict:
                tmp.append(sumDict[cID].copy()) # add each df to list
        df = pd.concat(tmp, ignore_index=True) # create one df
        df.loc[:, 'Values'] = df.Values.astype(int)
        df = df.groupby('Values').sum() # aggregate by raster value
        df = df.reset_index()
        df.rename(columns={'Count':'Acc'}, inplace=True) # already converted to area in agg step at catchment level
        df = df.merge(sumDict[curID], on='Values', how='outer') # add individiual catchment stats back in
        df.rename(columns={'Count':'Cat'}, inplace=True)
        df.fillna(0, inplace=True) # catchments without class get 0
        return df

####################################################################################################
#############################       GENERATE NETWORK      ##########################################
####################################################################################################
def generateDirectLinks(PLUS, direction):
    """
    Method: generateDirectLinks()
    Purpose: To compute the complete upstream network for the specified catchment
    Params:  PLUS - dataframe of the routing database
             direction - upstream or downstream
    Returns: allUp - dictionary of direct links
    """
    print("Find direct links...")
    #Create dict of direct upstream connections -- faster than referencing table directly each time
    allUp = {}
    allLines = list(set(list(PLUS['FL_ComID']))) #get all unique IDs in relation table
    for s in allLines:
        if direction == 'downstream':
            fromN = list(PLUS[PLUS['FL_ComID'] == s]['cfromnode'])[0]
            if -100 == fromN: #Nan values converted to -100; ignore these
                x = []
            else:
                x = list(PLUS[PLUS['ctonode'] == fromN]['FL_ComID'])
        else:
            fromN = list(PLUS[PLUS['FL_ComID'] == s]['ctonode'])[0]
            if -100 == fromN: #Nan values converted to -100; ignore these
                x = []
            else:
                x = list(PLUS[PLUS['cfromnode'] == fromN]['FL_ComID'])
        allUp[s] = x.copy()
    return allUp

def generateNetwork(allUp, s, PLUS, direction):
    """
    Method: generateNetwork()
    Purpose: To compute the complete upstream network for the specified catchment
    Params: allUp - dictionary containing the direct upstream connections for each unique
                    FL_ComID in the routing database
            s - the current catchment ID
            PLUS - dataframe of the routing database
    Returns: toSum - list of complete upstream connections (catchments and flowlines)
    """
    toSum = [int(s)] 
    new = [int(s)]
    while len(new) > 0:
        new, x, allUp = iterCheck(new, allUp, PLUS, direction) 
        if len(x) > 0:
            toSum = list(set(toSum) - set(x))
            new = list(set(new) - set(x))
        new = list(set(new) - set(toSum)) #remove any items from new that are already stored in the network
        toSum = toSum + new
    return list(set(toSum))

def iterCheck(coms, allUp, PLUS, direction):
    """
    Method: iterCheck()
    Purpose: find direct upstream connections for all FL_ComIDs in a list and build a list
            of the unique IDs found. Called by generateNetwork
    Params: coms - list of FL_ComIDs to find direct upstream connections of
            allUp - dictionary containing the direct upstream connections for each unique
                    FL_ComID in the routing database
            PLUS - dataframe of the routing database
    Returns: new - list of unique upstream catchment IDs
            notFound - list of unique IDs that did not exist in allUp
    """
    for c in coms:
        if c not in dirLinks: # ensure unique ID is in dict - this shouldn't execute
            if direction == 'downstream':
                fromN = list(PLUS[PLUS['FL_ComID'] == c]['cfromnode'])[0]
                if -100 == fromN: #Nan values converted to -100; ignore these
                    x = []
                else:
                    x = list(PLUS[PLUS['ctonode'] == fromN]['FL_ComID'])
            else:
                fromN = list(PLUS[PLUS['FL_ComID'] == c]['ctonode'])[0]
                if -100 == fromN: #Nan values converted to -100; ignore these
                    x = []
                else:
                    x = list(PLUS[PLUS['cfromnode'] == fromN]['FL_ComID'])
            x = [float(y) for y in x]
            dirLinks[c] = x.copy()
        if c in dirLinks:
            if len(dirLinks[c]) > 0:
                new = new + dirLinks[c]
        else:
            print(c, ": not in dirLinks")
            notFound.append(c)
    return list(set(new)), notFound, dirLinks

####################################################################################################
#########################       WRITE DATA       ###################################################
####################################################################################################
def writeShapefile(NHDPath, path, name, data, catchID):
    """
    Method: writeShapefile()
    Purpose: Produce the catchment statistics and accumulated statistics as a shapefile. Joins
            the statistics with the catchments.
    Params: NHDPath - path to the catchments
            path - directory to write the CSV to
            name - desired filename
            data - dataframe of statistics
            catchID - name of catchment ID field to merge on
    Returns: None
    """
    NHD = gpd.read_file(NHDPath) #read shapefile as geopandas df
    NHD = NHD[[catchID, 'geometry']] #was featureid
    #merge data
    NHD = NHD.merge(data, on=[catchID], how='left') #right instead of left to exclude writing catchments with no data
    f = os.path.join(path, name+'.shp')
    NHD.to_file(f)

def createTable(sumDict, rasVals):
    """
    Method: writeTable()
    Purpose: Produce the catchment statistics and accumulated statistics as a CSV file
    Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
            path - directory to write the CSV to
            name - desired filename (must end with .csv)
            rasVals - list of unique values for classed data (used as column names)
    Returns: None
    """
    if len(rasVals) == 0:
        rasVals = ['Stat', 'Area']
    cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
    for idx, v in enumerate(rasVals):
        cols[idx] = str(rasVals[idx])
        cols[idx+len(rasVals)] = "Ws_"+str(rasVals[idx])
    df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
    df['FEATUREID'] = df.index
    df = df.sort_values(by=['FEATUREID'])
    finCols = ['FEATUREID'] + cols
    df = df[finCols]
    df = df.fillna(0)
    return df

"""
Global variables needing user input
rasType: 0 (continuous raster) or 1 (classified raster) or 2 (tabular)
statType: "MIN", "MAX", "MEDIAN", "MEAN", or "SUM" (choose one statistic for continuous raster only)
convFact: conversion factor to apply to accumulated values; recommend output be in sq km
VPUzone: NHD zone that the catchments are in (02 is MidAtlantic)
catchID: Name of column with unique catchment ID in the catchments shapefile
tabID: if using tabular data containing catchment statistics, this is the name of the column with the catchment IDs
tabCols: if using tabular data containing catchment statistics, this is a list of columns to accumulate
"""
rasType = 0
statType = 'SUM'
catchID = 'FEATUREID'
outputCSV = "P6LU_v1_acres"
outputSHP = "P6LU_v1_acres"
OUTPUT_PATH = r"/media/imagery/1m_LU_2017/version1/NHD_Summaries"
convFact =  1 / 4046.83 # square meters to acres
direction = 'downstream' # 'upstream'
VPUzone = "02"
tabID = 'COMID'
tabCols = ['IR', 'INR', 'TCI', 'TG', 'TCT', 'FORE', 'WLF', 'WLO', 'WLT', 'MO', 'CRP', 'PAS', 'WAT']

#Data Paths
NHD_SHP = "/media/imagery/usgs_sc/smcdonald/Data/MD_healthy_h20sheds_upstream_COMIDs_albers.shp"
sassPath = r"/data/gis2/sm/NHD_Test/Version3_MWDB/Output/SassDB_02.csv"
# sassPath = r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/Test_Input/nhdplusv2_us.sas7bdat"
input_data_dict = { #input rasters or tabular catchment data
    # 'IMP17'   : r"/media/imagery/1m_LU_2017/version1/LUMM_10m/tot_imp_2017_v1_10m.tif",
    # 'NAT17'   : r"/media/imagery/1m_LU_2017/version1/LUMM_10m/tot_nat_2017_v1_10m.tif",
    'IR'    : r'/media/imagery/1m_LU_2017/version1/phase6_10m/IR_2017_10m.tif',
    'INR'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/INR_2017_10m.tif',
    'TCI'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/TCI_2017_10m.tif',
    'TG'    : r'/media/imagery/1m_LU_2017/version1/phase6_10m/TG_2017_10m.tif',
    'TCT'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/TCT_2017_10m.tif',
    'FORE'  : r'/media/imagery/1m_LU_2017/version1/phase6_10m/FOR_2017_10m.tif',
    'WLT'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/WLT_2017_10m.tif',
    'WLO'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/WLO_2017_10m.tif',
    'WLF'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/WLF_2017_10m.tif',
    'MO'    : r'/media/imagery/1m_LU_2017/version1/phase6_10m/MO_2017_10m.tif',
    'CRP'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/CRP_2017_10m.tif',
    'PAS'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/PAS_2017_10m.tif',
    'WAT'   : r'/media/imagery/1m_LU_2017/version1/phase6_10m/WAT_2017_10m.tif',
}

"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    invalidData = False #Flag to clean exit if any function fails
    #Read input data
    startTime = float(time.time())
    sumDict = {}
    allUp = {}
    PLUS = readRouting(sassPath, VPUzone)
    if len(PLUS) == 0:
        invalidData = True
    #Generate direct upstream network
    if not invalidData:
        allUp = generateDirectLinks(PLUS, direction)
        if len(allUp) == 0:
            invalidData = True   
    data_dict = {}
    if rasType == 0 or rasType == 1: #Raster data - continuous or categorical
        shapes, crs = readCatchmentGeoms(NHD_SHP, catchID)
        # verify all shapes keys are in PLUS
        p_ids = set(list(PLUS['FL_ComID']))
        s_ids = set(list(shapes.keys()))
        missing = s_ids - p_ids
        if len(missing) > 0:
            print("Geometry IDs do not exist in routing DB: ", len(missing))
            invalidData = True   
        if len(shapes) == 0:
            invalidData = True
        for ras in input_data_dict:
            print(ras)
            cellSize, nodata = readRas(input_data_dict[ras])
            if cellSize == 0:
                invalidData = True
            else:
                sumDict = srs.aggregateRas_mp(shapes, statType, input_data_dict[ras], rasType, convFact, cellSize, crs)
                if len(sumDict) == 0:
                    invalidData = True
                    print(f"{ras} sumDict empty")
                data_dict[ras] = sumDict.copy()
                if len(input_data_dict) > 1:
                    del sumDict
        if len(input_data_dict) > 1: # multiple rasters ran - convert and condense the data into 1 dict
            print("Condensing raster info into one dictionary...")
            df_list = []
            for ras in input_data_dict:
                df = pd.DataFrame.from_dict(data_dict[ras], orient='index', columns=[ras, 'TotAcres'])
                df['FEATUREID'] = df.index
                df = df[['FEATUREID', ras]]
                df_list.append(df.copy())
                del df
            df = pd.DataFrame()
            for d in df_list:
                if len(df) == 0:
                    df = d.copy()
                else:
                    df = df.merge(d, on='FEATUREID')
            del df_list
            df = df.reset_index()
            df = df.set_index('FEATUREID')
            df = df[list(input_data_dict.keys())]
            sumDict = df.to_dict(orient='index') 
            df.to_csv(f"{OUTPUT_PATH}/summary_stats.csv")
            del df
            for x in sumDict: # dict of dicts - convert to dict of lists
                sumDict[x] = list(sumDict[x].values())
            rasType = 1
            rasVals = list(input_data_dict.keys())
        sumDict = accumulate(rasType, rasVals, allUp, sumDict, convFact)
    elif rasType == 2: #tabular data
        sumDict = readTabular(inputData, tabID, tabCols, PLUS)
        # verify all shapes keys are in PLUS
        p_ids = set(list(PLUS['FL_ComID']))
        s_ids = set(list(sumDict.keys()))
        missing = s_ids - p_ids
        if len(missing) > 0:
            print("Geometry IDs do not exist in routing DB: ")
            for m in missing:
                print("\t", m)
            invalidData = True 
        if len(sumDict) == 0:
            invalidData = True
        else:
            sumDict = accumulate(rasType, tabCols, allUp, sumDict, convFact)
            if len(sumDict) == 0:
                invalidData = True

    lastTime = float(time.time())
    print("\nTime: ", str((lastTime- startTime)/60), " min\n")

    #write out data
    lastTime = float(time.time())
    if rasType == 2:
        rasVals = tabCols.copy()
    if not invalidData:
        data = createTable(sumDict, rasVals)
        if len(outputCSV) > 0:
            csv_path = os.path.join(OUTPUT_PATH, outputCSV+'.csv')
            data.to_csv(csv_path, index=False)
        if len(outputSHP) > 0:
            writeShapefile(NHD_SHP, OUTPUT_PATH, outputSHP, data, catchID)
    print("\nTotal Time: ", str((lastTime - startTime)/60), " min\n")
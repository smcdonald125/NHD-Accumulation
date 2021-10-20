"""
Purpose: Calculate raster statistics at the catchment level and accumulate the data. The input raster can be continuous data,
        where the following statistics can be calculated: MAX, MIN, MEAN, MEDIAN, SUM. The raster can also be classed data, where the
        total area per class is calculated. The catchment statistics and accumulated data can be written as a table to a CSV and
        can be written as a shapefile, joined to the catchments.
Data:   ecoSheds high res catchments as shapefile - Shapefile to act as zones; Is also used to build upstream network
        Raster File (tif) - user must specify if the data is continuous or classed
Note: The user is responsible for ensuring that the catchments and the raster are in the same projection.
Authors:  Sarah McDonald, Geographer, USGS
Contact: smcdonald@chesapeakebay.net or smcdonald@usgs.gov
"""

#Import Libs
import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.mask
import fiona
import math
import os
import time
import sys
import multiprocessing as mp
from shapely.geometry import Polygon

##########################################################################################################
#################################### READ DATA ###########################################################
##########################################################################################################
def readDataEco(catchPath, catchName):
    """
    Method: readDataEco()
    Purpose: Read ecosheds catchments and create "PlusFlow table" (table with catchment connections) and build
            dictionary of fiona geometries for each catchment.
    Params: catchPath - path to ecoSHEDs catchments
            catchName - string of catchment field ID (FEATUREID)
    Returns: PLUS - pandas dataframe of catchment connections
            shapes - dictionary of catchments geometries
    """
    shapes = {}
    PLUS_list = []
    if os.path.isfile(catchPath):
        PLUS = gpd.read_file(catchPath)[['FEATUREID', 'NextDownID']]
        with fiona.open(catchPath, "r") as geoms:
            for feature in geoms:
                shapes[feature['properties'][catchName]] = feature['geometry']  #dict of geomtries with FeatureID as key
    else:
        print("Catchment path does not exist: ", catchPath)
        sys.exit()

    return PLUS, shapes

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

##########################################################################################################
############################## ACCUMULATE DATA ###########################################################
##########################################################################################################
def generateDirectLinks(PLUS, toFromNames, direction):
    """
    Method: generateDirectLinks()
    Purpose: To compute the complete upstream or downstream network for the specified catchment
    Params:  PLUS - dataframe of the routing database
             toFromNames - list with column names; index 0 is downstream connection name, index 1 is current catchment name
             direction - upstream or downstream
    Returns: dirLinks - dict of direct network connections in the specified direction
    """
    #Create dict of direct connections -- faster than referencing table directly each time
    dirLinks = {}
    allLines = list(set(list(PLUS[toFromNames[0]] + PLUS[toFromNames[1]])))
    allLines = [float(x) for x in allLines]
    if 0 in allLines:
        allLines.remove(0)
    for cID in allLines:
        #find direct connections 
        if direction == 'downstream':   # accumulate with the flow direction
            x = list(PLUS[PLUS[toFromNames[0]] == float(cID)][toFromNames[1]]) # list of IDs that are upstream of current cID
        else: # accumulate against the flow direction
            x = list(PLUS[PLUS[toFromNames[1]] == float(cID)][toFromNames[0]]) # list of IDs that are downstream of current cID
        x = [float(i) for i in x]
        x = list(set(x))
        if 0 in x:
            x.remove(0)
        dirLinks[float(cID)] = x.copy()
    return dirLinks

def generateNetwork(dirLinks, s, PLUS, direction, toFromNames):
    """
    Method: generateNetwork()
    Purpose: To compute the complete up/down stream network for the specified catchment
    Params: dirLinks - dictionary containing the direct network connections in the specified direction for each unique
                    ID in the NHD PlusFLow / ecsheds catchments or flowlines dbf
            s - the current catchment ID
            PLUS - dataframe of all connections
            direction - string, either 'upstream' or 'downstream'
            toFromNames - list with column names; index 0 is downstream connection name, index 1 is current catchment name
    Returns: toSum - list of complete up/down stream connections (catchments and flowlines)
    """
    toSum = [float(s)]
    new = [float(s)]
    while len(new) > 0:
        new, x, dirLinks = iterCheck(new, dirLinks, PLUS, direction, toFromNames) #iterative version of recursive check()
        if len(x) > 0:
            toSum = list(set(toSum) - set(x))
            new = list(set(new) - set(x))
        new = list(set(new) - set(toSum)) #remove any items from new that are already stored in the network
        toSum = toSum + new
    return list(set(toSum))

def iterCheck(coms, dirLinks, PLUS, direction, toFromNames):
    """
    Method: iterCheck()
    Purpose: find direct connections for all IDs in a list and build a list
            of the unique IDs found. Called by generateNetwork
    Params: coms - list of IDs to find direct connections of
            dirLinks - dictionary containing the direct connections in the sprecified direction for each unique
                    ID in the NHD PlusFLow / ecsheds catchments or flowlines dbf
            PLUS - dataframe of all connections
            direction - string of upstream or downstream
            toFromNames - list with column names; index 0 is downstream connection name, index 1 is current catchment name
    Returns: new - list of unique upstream catchment IDs
            notFound - list of unique IDs that did not exist in dirLinks
            dirLinks - original dictionary, can be updated in this method so need to return any changes
    """
    new = []
    notFound = []
    for c in coms:
        if c not in dirLinks: # ensure unique ID is in dict - this shouldn't execute
            if direction == 'downstream':
                x = list(PLUS[PLUS[toFromNames[0]] == c][toFromNames[1]])
            else:
                x = list(PLUS[PLUS[toFromNames[1]] == c][toFromNames[0]])
            x = [float(y) for y in x]
            dirLinks[c] = x.copy()
        if c in dirLinks:
            if len(dirLinks[c]) > 0:
                new = new + dirLinks[c]
        else:
            print(c, ": not in dirLinks")
            notFound.append(c)
    return list(set(new)), notFound, dirLinks

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
    allCatch, PLUS, dirLinks, network, sumDict, convFact, direction, rasType, toFromNames = args
    #Accumulate statistic for each catchment
    allSum = list(sumDict)
    newSum, newUp = {}, {}
    for s in allCatch:
        if s not in network:
            network[s] = generateNetwork(dirLinks, s, PLUS, direction, toFromNames) #all network connections
            newUp[s] = network[s].copy()
        upCatch = list((set(network[s])&set(allSum))) #all network catchments that are in sumDict
        newSum[s] = accumulate_values(sumDict, upCatch, convFact, s, rasType) #list of catchment values (original scale) and accumulated values (converted scale)
    return newSum, newUp #return all data in dict

def accumulate_mp(rasType, PLUS, dirLinks, network, sumDict, convFact, direction, toFromNames):
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
            toFromNames - list with column names; index 0 is downstream connection name, index 1 is current catchment name
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
            network - dictionary storing full network in specified direction
    """
    cpus_minus_1 = mp.cpu_count() - 2
    if cpus_minus_1 == 0:
        cpus_minus_1 = 1

    allCatch = list(sumDict)
    batch_size = len(allCatch) / cpus_minus_1
    if batch_size % 1 != 0:
        batch_size = int((batch_size) + 1)
    else:
        batch_size = int(batch_size)
    
    chunk_iterator = []
    for i in range(cpus_minus_1):
        mn, mx = i * batch_size, (i + 1) * batch_size
        keys = allCatch[mn:mx]
        gdf_args = keys, PLUS, dirLinks, network, sumDict, convFact, direction, rasType, toFromNames
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

##########################################################################################################
############################## SUMMARIZE RASTER ##########################################################
##########################################################################################################
def aggregateRas(args):
    """
    Method: aggregateRas()
    Purpose: calculate raster statistics for each catchment, find the upstream network for each catchment and
            accumulate the raster statistics.
    Params: keys - list of catchment IDs
            geoms - list of catchment geometries
            rasPath - path to the raster data
            rasType - 0 or 1 specifying if raster is continuous or classed
            stat - statistic type to be calculated for the catchments (only needed for continuous data)
            convFact - user defined conversion factor to apply to raster area
            cellSize - raster cell sized; needed to convert pixel count to area for classed data
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
    """
    keys, geoms, rasPath, rasType, stat, convFact, cellSize = args
    sumDict = {}
    #Calculate statistic for each catchment
    with rio.open(rasPath) as src:
        noData = src.nodatavals[0]
        if rasType == 0: #continuous raster
            for s in range(len(keys)):
                val, ar = agg_cont(src, geoms[s], stat, convFact, noData)
                if ar != -1:
                    sumDict[keys[s]] = [float(val), float(ar)]
        elif rasType == 1: #classed raster
            for s in range(len(keys)):
                df = agg_classed(src, geoms[s], noData, cellSize, convFact)
                if len(df) > 0: # replace with .empty?
                    sumDict[keys[s]] = df.copy()
                    del df
    return sumDict 

def aggregateRas_mp(shapes, stat, rasPath, rasType, convFact):
    """
    Method: aggregateRas_mp()
    Purpose: Divide catchments up to be multi-processed through statistic calculations.
    Params: shapes - dict of catchment IDs and geometries
            stat - statistic type to be calculated for the catchments (only sued for continuous data)
            rasPath - path to the raster data
            rasType - 0 or 1 specifying if raster is continuous or classed
            convFact - use defined conversion factor to apply to raster area
            cellSize - raster cell sized; needed to convert pixel count to area for classed data
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
    """
    rt = time.time()
    cpus_minus_1 = mp.cpu_count() - 2
    if cpus_minus_1 == 0:
        cpus_minus_1 = 1

    pool = mp.Pool(processes=cpus_minus_1)

    batch_size = len(shapes) / cpus_minus_1
    if batch_size % 1 != 0:
        batch_size = int((batch_size) + 1)
    else:
        batch_size = int(batch_size)

    chunk_iterator = []
    for i in range(cpus_minus_1): # was num_chunks
        mn, mx = i * batch_size, (i + 1) * batch_size
        keys = list(shapes)[mn:mx]
        geoms = list(shapes.values())[mn:mx]
        gdf_args = keys, geoms, rasPath, rasType, stat, convFact, cellSize
        chunk_iterator.append(gdf_args)

    results = pool.map(aggregateRas, chunk_iterator) #list of dictionaries
    pool.close()

    result = {} #create one dict from list of dicts
    for d in results:
        result.update(d)
    del results

    return result
         
def agg_classed(src, geom, noData, cellSize, convFact):
    """
    Method: agg_classed()
    Purpose: Mask the classed raster to the current catchment and calculate the cell count
            for each class. Store data in a dataframe.
    Params: src - Rasterio open DatasetReader object for the raster
            geom - fiona geometry of the catchment
            noData - raster noData value
            cellSize - raster cell sized; needed to convert pixel count to area for classed data
            convFact - user defined conversion factor to apply to area
    Returns: df - pandas dataframe of raster values and counts
    """
    try:
        ary, t = rio.mask.mask(src, [geom], crop=True, all_touched=False) #mask by current catchment
        vals, counts = np.unique(ary, return_counts=True) 
        vals = list(vals)
        convFact = convFact * cellSize * cellSize # conversion factor times cell area
        counts = [x * convFact for x in list(counts)] # convert count to area
        if noData in vals:
            i = vals.index(noData)
            vals.remove(noData) #update to ignore no data value -- 255 hard-coded for NLCD
            del counts[i]
        if len(vals) > 0:
            df = pd.DataFrame(data={'Values':vals, 'Count': counts})
            return df
        else:
            return pd.DataFrame()
    except:
        return pd.DataFrame()

def agg_cont(src, geom, stat, convFact, noData):
    """
    Method: agg_cont()
    Purpose: Mask the continuous raster to the current catchment and calculate the specified statistic
            and raster area within the catchment. The statistic options are: MAX, MIN, MEAN, MEDIAN and SUM.
    Params: src - Rasterio open DatasetReader object for the raster
            geom - fiona geometry of the catchment
            stat - string denoting the statistic type
            convFact - conversion factor to be applied to the area
            noData - raster noData value 
    Returns: catchSum - catchment statistic
            area - catchment area
    """
    try:
        ary, t = rio.mask.mask(src, [geom], crop=True, all_touched=False) #mask by current catchment
        area = ((ary != noData).sum())* convFact # number of pixels that are not no data
        #set no data value to nan for float rasters and 0 for integer rasters
        if type(ary[0][0][0]) == np.float32 or type(ary[0][0][0]) == np.float64:
            ary[ary == noData] = np.nan #set nodata to nan 
        else:
            ary[ary == noData] = 0
        if statType == "MAX":
            catchSum = float(np.amax(ary)) 
        elif statType == "MIN":
            catchSum = float(np.amin(ary)) 
        elif statType == "MEAN":
            catchSum = float(np.nanmean(ary)) 
        elif statType == "MEDIAN":
            catchSum = float(np.nanmedian(ary)) 
        elif statType == "SUM":
            catchSum = float(np.sum(ary))
        return catchSum, area #returns value
    except:
        return -1, -1
    return 0

##########################################################################################################
############################## WRITE DATA ################################################################
##########################################################################################################
def createTable(sumDict, catchID, rasType, statType):
    """
    Method: createTable()
    Purpose: Produce the catchment statistics and accumulated statistics as a dataframe file
    Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
            path - directory to write the CSV to
            name - desired filename (must end with .csv)
            catchID - name of catchment ID field
            rasType - 0 or 1 specifying if raster is continuous or classed
            statType - string stat for continuous data (SUM, MEDIAN, etc.)
    Returns: df - dataframe of data
    """
    if rasType == 0: #Continuous data has 1 stat and area
        rasVals = [statType, 'Area']
        cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
        for idx, v in enumerate(rasVals):
            cols[idx] = str(rasVals[idx])
            cols[idx+len(rasVals)] = "Ws_"+str(rasVals[idx])
        df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
        df[catchID] = df.index
        df = df.sort_values(by=[catchID])
        finCols = [catchID] + cols
        df = df[finCols]
    else: # Classed data
        allData = []
        for s in sumDict:
            tmp = sumDict[s].copy() # read in df
            vals = list(tmp['Values'])
            cols = [0] * (len(vals) * 2)
            for idx, v in enumerate(vals):
                cols[idx] = str(vals[idx])
                cols[idx+len(vals)] = "Ws_"+str(vals[idx])
            df = pd.DataFrame(data={'cols':[catchID]+cols}) # put columns in a single column for now
            df.loc[:, 'count'] = [s] + list(tmp['Cat']) + list(tmp['Acc'])
            df = df.set_index('cols').T #columns are now columns, ID and values are completed
            allData.append(df.copy())
        df = pd.concat(allData, ignore_index=True)
        df = df.sort_values(by=[catchID])
        df.fillna(0, inplace=True)
        cols = [x for x in list(df) if x != catchID]
        cols.sort() # organize columns
        cols = [catchID] + cols
        df = df[cols]
    return df

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
    NHD = NHD.merge(data, on=[catchID], how='right') #right instead of left to exclude writing catchments with no data
    f = os.path.join(path, name+'.shp')
    NHD.to_file(f)

##########################################################################################################
################################### HELPERS ##############################################################
##########################################################################################################
def validateInput(rasType, statType, outputCSV, outputSHP, OUTPUT_PATH, convFact, direction, NHD_SHP, RAS):
    """
    Method: validateInput()
    Purpose: Ensure user input is valid before running accumulation.
    """
    valid, outFlag, writeFlag = True, True, False
    if rasType not in [0,1]:
        print("INPUT ERROR: rasType must be either 0 or 1\n\t0: Continuous\n\t1: Categorical")
        valid = False
    if rasType == 0 and statType not in ["MIN", "MAX", "MEDIAN", "MEAN", "SUM"]:
        print("INPUT ERROR: statType for a continuous raster must be\n\tMIN, MAX, MEDIAN, MEAN, or SUM")
        valid = False
    if not os.path.isdir(OUTPUT_PATH):
        try:
            os.mkdir(OUTPUT_PATH)
            print(f"Created Output Directory: {OUTPUT_PATH}")
        except:
            print(f"INPUT ERROR: OUTPUT_PATH directory does not exist and failed to be created\n\t{OUTPUT_PATH}")
            valid = False
            outFlag = False
    if len(outputCSV) > 0:
        writeFlag = True
        if '.' in outputCSV: # remove .csv
            print(f"INPUT ERROR: outputCSV file includes extension\n\t{outputCSV}")
            valid = False
        elif outFlag:
            if os.path.isfile(os.path.join(OUTPUT_PATH, outputCSV+'.csv')):
                print(f"WARNING: output csv file already exists and will be overwritten\n\t{os.path.join(OUTPUT_PATH, outputCSV+'.csv')}")
    if len(outputSHP) > 0:
        writeFlag = True
        if '.shp' in outputSHP: # remove .csv
            print(f"INPUT ERROR: outputSHP file includes extension\n\t{outputSHP}")
            valid = False
        elif outFlag:
            if os.path.isfile(os.path.join(OUTPUT_PATH, outputSHP+'.shp')):
                print(f"WARNING: output shapefile already exists and will be overwritten\n\t{os.path.join(OUTPUT_PATH, outputSHP+'.shp')}")
    if not writeFlag:
        print("INPUT ERROR: outputCSV and outputSHP are blank\n\tUser must specificy output file name for at least one")
        valid = False
    if convFact <= 0:
        print("INPUT ERROR: convFact must be greater than 0\n\tNote: enter 1 for no conversion")
        valid = False
    if direction not in ['upstream', 'downstream']:
        print("INPUT ERROR: direction must be upstream or downstream")
        valid = False
    if not os.path.isfile(NHD_SHP):
        print(f"INPUT ERROR: NHD_SHP does not exist\n\t{NHD_SHP}")
        valid = False
    if not os.path.isfile(RAS):
        print(f"INPUT ERROR: RAS does not exist\n\t{RAS}")
        valid = False
    if not valid:
        sys.exit(1)

def printTime(message, prevTime):
    """
    Method: printTime()
    Purpose: Accept previous time and calculate difference with current time. Print message and time difference
            with approopriate time metric (sec, min, or hours).
    Params: message - string containing message to print before time elapsed
            prevTime - float of time.time() of time before the running the code you want to time
    Returns: curTime - current time used for time elapsed calculation
    """
    curTime = float(time.time())
    timeDif = float(curTime - prevTime)
    t = ' sec'
    if timeDif > 60:
        timeDif = timeDif / 60
        t = ' min'
        if timeDif > 60:
            timeDif = timeDif / 60
            t = ' hours'
    print(message, str(timeDif), t, "\n")
    return curTime
    

"""
Global variables needing user input
rasType: 0 (continuous raster) or 1 (categorical raster)
statType: "MIN", "MAX", "MEDIAN", "MEAN", or "SUM" (choose one statistic for continuous raster only)
outputCSV: string name for accumulated output in CSV format - do not include .csv extension
outputSHP: string name for accumulated output in shapefile format - do not include .shp extension
OUTPUT_PATH: string path to folder to write output to
convFact: conversion factor to apply to accumulated values; recommend output be in sq km or hectares
catchIDName: name of the columns denoting the unique catchment IDs
toFromPlusNames: list of length to denoting the column names for the to/from columns;
                 for NHD high res data ['ToNHDPID', 'FromNHDPID']
                 for ecoSHEDs ['NextDownID', 'FEATUREID'] 
"""
rasType = 0
statType = "SUM"
outputCSV = "Potomac_Eco_INR17_10m_upstream_20211019"
outputSHP = "Potomac_Eco_INR17_10m_upstream_20211019"
OUTPUT_PATH = r"G:/ImageryServer/usgs_sc/smcdonald/NHDv2_Agg/24k_Beta/ecoSheds/Potomac_20211019/" #/media/imagery
convFact =  1 / 4046.83 # square meters to acres
direction = 'upstream' # 'upstream'
#Data Paths
NHD_SHP = r"G:/ImageryServer/usgs_sc/smcdonald/NHDv2_Agg/24k_Beta/ecoSheds/Catchments_POT.shp" #/media/imagery  /spatial_02/Catchments02.shp
RAS = r"G:/ImageryServer/1m_LU_2017/version1/phase6_10m/INR_2017_10m.tif"

# These will remain constant for ecoSHEDs data
catchIDName = 'FEATUREID'
toFromPlusNames = ['NextDownID', 'FEATUREID'] 

"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    # Check user input
    validateInput(rasType, statType, outputCSV, outputSHP, OUTPUT_PATH, convFact, direction, NHD_SHP, RAS)
    print("Validated Input")
    #Read input data
    startTime = float(time.time())
    sumDict = {}
    network = {} #build once - code loop to run multiple datasets
    PLUS, shapes = readDataEco(NHD_SHP, catchIDName)
    prevTime = printTime('Read EcoSHEDs Time', startTime)
    invalidData = False
    if len(PLUS) == 0:
        invalidData = True
        print("No records in PLUS - exiting")
    if len(shapes) == 0:
        invalidData = True
        print("No catchment geometries in dictionary - exiting")
    #Generate network
    if not invalidData:
        dirLinks = generateDirectLinks(PLUS, toFromPlusNames, direction)
        prevTime = printTime('Generate Network Time', prevTime)
        if len(dirLinks) == 0:
            invalidData = True
            print("Generate Network failed - no records in dirLinks - exiting")
    if not invalidData:     
        cellSize, noData = readRas(RAS)
        print("Cell Size: ", cellSize)
        print("No Data: ", noData)
        prevTime = printTime('Read Raster Information Time', prevTime)
        if cellSize == 0:
            invalidData = True
            print("Cell size is 0 - exiting")
        else:
            sumDict = aggregateRas_mp(shapes, statType, RAS, rasType, convFact)
            prevTime = printTime('Summarize Raster Information Time', prevTime)
            if len(sumDict) == 0:
                invalidData = True
                print("No Catchment Summaries in sumDict - exiting")
            else:
                convFact = convFact * cellSize * cellSize
                print("Conversion Factor: ", convFact)
                sumDict, network = accumulate_mp(rasType, PLUS, dirLinks, network, sumDict, convFact, direction, toFromPlusNames)
                prevTime = printTime('Accumulate Time', prevTime)

    #write out data
    if not invalidData:
        data = createTable(sumDict, catchIDName, rasType, statType)
        if len(outputCSV) > 0:
            csv_path = os.path.join(OUTPUT_PATH, outputCSV+'.csv')
            data.to_csv(csv_path, index=False)
        if len(outputSHP) > 0:
            writeShapefile(NHD_SHP, OUTPUT_PATH, outputSHP, data, catchIDName)
    prevTime = printTime('Write Results Time', prevTime)

    printTime('Total Time', startTime)
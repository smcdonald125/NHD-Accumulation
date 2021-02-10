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


##########################################################################################################
#################################### READ DATA ###########################################################
##########################################################################################################

"""
Method: readDataEco()
Purpose: Read ecosheds catchments and create "PlusFlow table" (table with catchment connections) and build
         dictionary of fiona geometries for each catchment.
Params: catchPath - path to ecoSHEDs catchments
        catchName - string of catchment field ID (FEATUREID)
Returns: PLUS - pandas dataframe of catchment connections
         shapes - dictionary of catchments geometries
"""
def readDataEco(catchPath, catchName):
    shapes = {}
    if os.path.isfile(catchPath):
        PLUS = gpd.read_file(catchPath)
        PLUS = PLUS[['FEATUREID', 'NextDownID']]
        
        with fiona.open(catchPath, "r") as geoms:
            shapes = {feature['properties'][catchName]:feature['geometry'] for feature in geoms} #dict of geomtries with FeatureID as key
    else:
        print("Catchment path does not exist: ", catchPath)
        sys.exit()

    return PLUS, shapes


"""
Method: readRas()
Purpose: Get cell size and unique class values from categorical raster
Params:  rasPath - path to raster
         rasType - 0 (continuous raster) or 1 (classified raster)
Returns: cellSize - cell size of the raster
         rasVals - list of unique raster values for classed data (empty list for continuous)
"""
def readRas(rasPath, rasType):
    print("Reading raster data...")
    if os.path.isfile(rasPath):
        with rio.open(rasPath, 'r') as src:
            noData = src.nodatavals
            noData = noData[0]
            t = src.transform
            cellSize = t.a
            rasVals = []
            if rasType == 1: #if categorical data
                src_array = src.read(1)
                rasVals = list(np.unique(src_array))
                if noData in rasVals:
                    rasVals.remove(noData)

        return cellSize, rasVals, noData     
    else:
        print("Invalid file path: ", rasPath)
        sys.exit()


##########################################################################################################
############################## ACCUMULATE DATA ###########################################################
##########################################################################################################

"""
Method: generateUpstream()
Purpose: To compute the complete upstream network for the specified catchment
Params:  PLUS - dataframe of the routing database
Returns: toSum - list of complete upstream connections (catchments and flowlines)
"""
def generateUpstream(PLUS, toFromNames):
    print("Finding direct upstream links...")
    #Create dict of direct upstream connections -- faster than referencing table directly each time
    allUp = {}
    allLines = list(set(list(PLUS[toFromNames[0]] + PLUS[toFromNames[1]])))
    allLines = [float(x) for x in allLines]
    if 0 in allLines:
        allLines.remove(0)
    for cID in allLines:
        #find direct upstream connections    
        x = list(PLUS[PLUS[toFromNames[0]] == float(cID)][toFromNames[1]])
        x = [float(i) for i in x]
        x = list(set(x))
        if 0 in x:
            x.remove(0)
        allUp[float(cID)] = x
    return allUp


"""
Method: accumulate()
Purpose: accumulate the raster statistics
Params: convFact - the conversion factor specified by the user
        rasType - 0 or 1 specifying if raster is continuous or classed
        rasVals - list containing unique class values for a classed raster (empty list for continuous)
        allUp - dictionary of direct upstream links
        sumDict - dictionary where the catchment ID is the key and the data is a list of catchment statistics
Returns: sumDict - dictionary where the catchment ID is the key and the data is a list of catchment statistics and
                   accumulated statistics
"""
def accumulate(rasType, rasVals, allUp, upstream, sumDict, convFact): #switched allUp with upstream (already built full lists
#Accumulate statistic for each catchment
    print("Accumulating...")
    allCatch = list(sumDict.keys())
    if rasType == 0: #if continuous data -- update list of "unique values" to be the statistic and area
        rasVals = ['Stat', 'Area']
    numClasses = len(rasVals)
    for s in allCatch:
        if s not in upstream:
            upstream[s] = allUpstream(allUp, s, PLUS) #all upstream connections
        upCatch = list((set(upstream[s])&set(allCatch))) #all upstream catchments
        sumDict[s] = getSumList(sumDict, upCatch, numClasses, convFact, s) #list of catchment values (original scale) and accumulated values (converted scale)
    return sumDict, upstream #return all data in dict


"""
Method: getSumList()
Purpose: Sum each statistic value for all catchments in the upstream network.
Params: sumDict - dictionary containing catchment statistics
        toSum - list of upstream catchment IDs
        numClasses - integer denoting the number of statistics to sum (# of classes for classed data
                    or 2 for continuous[statistic, area])
        convFact - the conversion factor to be applied to the accumulated values 
        curID - the "current" catchment ID
Returns: allStats - list of catchment statistics and accumulated statistics for the current catchment
"""
def getSumList(sumDict, toSum, numClasses, convFact, curID):
    sums = [0] * numClasses
    for cID in toSum:
        if cID in sumDict: #check that upstream catchment is in the dict
            tmp = sumDict[cID] #get list of all stats for upstream catchment
            for idx in range(numClasses): #for each class -- sum the data
                sums[idx] += tmp[idx]
    for i in range(len(sums)):
        sums[i] = sums[i] * convFact #convert data
    allStats = sumDict[curID] + sums #create full list in the order "catchment statistics for each class and accumulated stats for each class"
    return allStats

"""
Method: allUpstream()
Purpose: To compute the complete upstream network for the specified catchment
Params: allUp - dictionary containing the direct upstream connections for each unique
                ID in the NHD PlusFLow / ecsheds catchments or flowlines dbf
        s - the current catchment ID
        PLUS - dataframe of the connections
        
Returns: toSum - list of complete upstream connections (catchments and flowlines)
"""
def allUpstream(allUp, s, PLUS):
    toSum = [float(s)] #was allUp[s] + [s] changed since moving allUp to iterCheck
    new = [float(s)]
    while len(new) > 0:
        new, x = iterCheck(new, allUp, PLUS) #iterative version of recursive check()
        if len(x) > 0:
            toSum = list(set(toSum) - set(x))
            new = list(set(new) - set(x))
        new = list(set(new) - set(toSum)) #remove any items from new that are already stored in the network
        toSum = toSum + new
    return list(set(toSum))

"""
Method: iterCheck()
Purpose: find direct upstream connections for all IDs in a list and build a list
         of the unique IDs found. Called by allUpstream
Params: coms - list of IDs to find direct upstream connections of
        allUp - dictionary containing the direct upstream connections for each unique
                ID in the NHD PlusFLow / ecsheds catchments or flowlines dbf
        PLUS - dataframe of the connections
Returns: new - list of unique upstream catchment IDs
         notFound - list of unique IDs that did not exist in allUp
"""
def iterCheck(coms, allUp, PLUS):
    new = []
    notFound = []
    for c in coms:
        if c not in allUp:
            x = list(PLUS[PLUS['NextDownID'] == c]['FEATUREID'])
            x = [float(y) for y in x]
            allUp[c] = x
        if c in allUp:
            if len(allUp[c]) > 0:
                new = new + allUp[c]
        else:
            print(c, ": not in allUp")
            notFound.append(c)
    return list(set(new)), notFound

##########################################################################################################
############################## SUMMARIZE RASTER ##########################################################
##########################################################################################################

"""
Method: aggregateRas()
Purpose: calculate raster statistics for each catchment, find the upstream network for each catchment and
         accumulate the raster statistics.
Params: shapes - dict of fiona geometries for the catchments
        stat - statistic type to be calculated for the catchments (only needed for continuous data)
        rasPath - path to the raster data
        cellSize - raster cell size
        rasType - 0 or 1 specifying if raster is continuous or classed
        rasVals - list containing unique class values for a classed raster (empty list for continuous)
Returns: sumDict - dictionary where the catchment ID is the key and the data is a list of catchment statistics and
                   accumulated statistics
"""
def aggregateRas(shapes, stat, rasPath, cellSize, rasType, rasVals, noData):
    print("Calculating Catchment Statistics...")
    sumDict = {}
    #Calculate statistic for each catchment
    with rio.open(rasPath) as src:
        if rasType == 0: #continuous raster
            for s in shapes:
                val, ar = getSum(src, shapes[s], stat, cellSize, noData)
                sumDict[s] = [float(val), float(ar)]
        elif rasType == 1: #classed raster
            for s in shapes:
                vals = getSumDis(src, shapes[s], noData, rasVals)
                sumDict[s] = vals
    return sumDict

"""
Method: getSumDis()
Purpose: Mask the classed raster to the current catchment and calculate the cell count
         for each class. Organize the counts in a list that will align with the total 
         number of unique classes found in the whole raster.
Params: src - Rasterio open DatasetReader object for the raster
        geom - fiona geometry of the catchment
        noData - raster noData value
        rasVals - list of all unique classes in the raster
Returns: finVals - list of class pixel counts within the catchment
"""         
def getSumDis(src, geom, noData, rasVals):
    try:
        ary, t = rio.mask.mask(src, [geom], crop=True, all_touched=False) #mask by current catchment
        vals, counts = np.unique(ary, return_counts=True) 
        vals = list(vals)
        counts = list(counts)
        if noData in vals:
            i = vals.index(noData)
            vals.remove(noData) #update to ignore no data value -- 255 hard-coded for NLCD
            del counts[i]
        finVals = [0 for x in range(len(rasVals))] #empty list with same length as unique raster values
        for idx, v in enumerate(vals): #for each class found within the catchment
            finVals[rasVals.index(v)] = counts[idx] #set count of unique class in same order as the unique raster values
        return  finVals #returns list of counts
    except:
        print("Could Not Mask Raster -- Check Projections of Files - exiting")
        sys.exit()

"""
Method: getSum()
Purpose: Mask the continuous raster to the current catchment and calculate the specified statistic
         and raster area within the catchment. The statistic options are: MAX, MIN, MEAN, MEDIAN and SUM.
Params: src - Rasterio open DatasetReader object for the raster
        geom - fiona geometry of the catchment
        statType - string denoting the statistic type
        convFact - conversion factor to be applied to the area
        noData - raster noData value 
Returns: catchSum - catchment statistic
         area - catchment area
"""
def getSum(src, geom, stat, convFact, noData):
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
        print("Could Not Mask Raster -- Check Projections of Files - exiting")
        sys.exit()
    return 0

##########################################################################################################
############################## WRITE DATA ################################################################
##########################################################################################################

"""
Method: writeTable()
Purpose: Produce the catchment statistics and accumulated statistics as a CSV file
Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
        path - directory to write the CSV to
        name - desired filename (must end with .csv)
        rasVals - list of unique values for classed data (used as column names)
Returns: None
"""
def writeTable(sumDict, path, name, rasVals, catchID):
    print("Writing CSV...")
    if name[-4:] == ".csv":
        if len(rasVals) == 0:
            rasVals = ['Stat', 'Area']
        cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
        for idx, v in enumerate(rasVals):
            cols[idx] = str(rasVals[idx])
            cols[idx+len(rasVals)] = "Ws_"+str(rasVals[idx])
        df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
        df[catchID] = df.index
        df = df.sort_values(by=[catchID])
        finCols = [catchID] + cols
        df = df[finCols]
        fil = os.path.join(path, name)
        df.to_csv(fil, index=False)
    else:
        print("Output File Must be of Type .csv", name)

"""
Method: writeShapefile()
Purpose: Produce the catchment statistics and accumulated statistics as a shapefile. Joins
         the statistics with the NHD+ catchments.
Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
        NHDPath - path to the catchments
        path - directory to write the CSV to
        name - desired filename (must end with .csv)
        rasVals - list of unique values for classed data (used as column names)
Returns: None
"""
def writeShapefile(sumDict, NHDPath, path, name, rasVals, catchID):
    # Read catchments as geodataframe -- used to join with final results
    if os.path.isfile(NHDPath):
        NHD = gpd.read_file(NHDPath) #read shapefile as geopandas df
        NHD = NHD[[catchID, 'geometry']] #was featureid
        print("Writing Shapefile...")
        if name[-4:] == ".shp":
            f = os.path.join(path, name)
            if len(rasVals) == 0: #continous data
                rasVals = ['Stat', 'Area']
            cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
            for idx, v in enumerate(rasVals):
                cols[idx] = str(rasVals[idx])
                cols[idx+len(rasVals)] = "Ws_"+str(rasVals[idx])
            df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
            df[catchID] = df.index
            NHD = NHD.merge(df, on=[catchID], how='left')
            NHD.to_file(f)
        else:
            print("Output File Must be of Type .shp", name)
    else:
        print("Invalid File Path: ", NHDPath)
        print("No catchment shapefile to join with data")


"""
Method: printTime()
Purpose: Accept previous time and calculate difference with current time. Print message and time difference
         with approopriate time metric (sec, min, or hours).
Params: message - string containing message to print before time elapsed
        prevTime - float of time.time() of time before the running the code you want to time
Returns: curTime - current time used for time elapsed calculation
"""
def printTime(message, prevTime):
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
outputType: [1, 1] (first index is output type CSV and second is shapefile; 0 is do NOT write and 1 is write)
convFact: conversion factor to apply to accumulated values; recommend output be in sq km or hectares
catchIDName: name of the columns denoting the unique catchment IDs
toFromPlusNames: list of length to denoting the column names for the to/from columns;
                 for NHD high res data ['ToNHDPID', 'FromNHDPID']
                 for ecoSHEDs ['NextDownID', 'FEATUREID'] 
"""
rasType = 0
statType = "SUM"
outputType = [1, 1]
convFact = 1e-6 # this value will be multiplied by cell size squared automatically later; this example converts m2 to km2
catchIDName = 'FEATUREID'
toFromPlusNames = ['NextDownID', 'FEATUREID'] 

"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    #Data Paths
    NHD_SHP = r"" #path to catchments shapefile
    RAS = r"" #path to raster (tif) to summarize
    OUTPUT_PATH = r"" #path to location to write accumulated data (do not include file name)

    #Read input data
    startTime = float(time.time())
    sumDict = {}
    upstream = {} #build once - code loop to run multiple datasets
    PLUS, shapes = readDataEco(NHD_SHP, catchIDName)
    prevTime = printTime('Read EcoSHEDs Time', startTime)
    invalidData = False
    if len(PLUS) == 0:
        invalidData = True
        print("No records in PLUS - exiting")
    #Generate upstream network
    if not invalidData:
        allUp = generateUpstream(PLUS, toFromPlusNames)
        prevTime = printTime('Generate Upstream Time', prevTime)
        if len(allUp) == 0:
            invalidData = True
            print("Generate Upstream failed - no records in allUp - exiting")
    if not invalidData:     
        cellSize, rasVals, noData = readRas(RAS, rasType)
        prevTime = printTime('Read Raster Information Time', prevTime)
        if cellSize == 0:
            invalidData = True
            print("Cell size is 0 - exiting")
        else:
            if len(shapes) == 0:
                invalidData = True
                print("No catchment geometries in dictionary - exiting")
            else:
                sumDict = aggregateRas(shapes, statType, RAS, cellSize, rasType, rasVals, noData)
                prevTime = printTime('Summarize Raster Information Time', prevTime)
                if len(sumDict) == 0:
                    invalidData = True
                    print("No Catchment Summaries in sumDict - exiting")
                else:
                    convFact = convFact * cellSize * cellSize
                    sumDict, upstream = accumulate(rasType, rasVals, allUp, upstream, sumDict, convFact) #return upstream - only create once - future update loop ras summaries
                    prevTime = printTime('Accumulate Time', prevTime)

    #write out data
    if not invalidData:
        if outputType[0] == 1:
            writeTable(sumDict, OUTPUT_PATH, "POT_IMP_ECO_021021.csv", rasVals, catchIDName)
        if outputType[1] == 1:
            writeShapefile(sumDict, NHD_SHP, OUTPUT_PATH, "POT_IMP_ECO_021021.shp", rasVals, catchIDName)

    prevTime = printTime('Write Results Time', prevTime)

    printTime('Total Time', startTime)

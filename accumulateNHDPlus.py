"""
Purpose:Calculate raster statistics at the NHD+ catchment level and accumulate the data. The input raster can be continuous data,
        where the following statistics can be calculated: MAX, MIN, MEAN, MEDIAN, SUM. The raster can also be classed data, where the
        total area per class is calculated. The catchment statistics and accumulated data can be written as a table to a CSV and
        can be written as a shapefile, joined to the NHD+ catchments. The upstream networks are built using a routing database, produced
        by Mike Wieczorek. The routing database can be reduced to the VPU zone level, and can be read in as its original form (SASS database)
        or a CSV file. 
Data:   NHD+ catchments as shapefile - Shapefile to act as zones - must have FEATUREID and AreaSqKM fields
        Routing Database as SASS database or CSV- Must have FL_ComID, cfromnode, ctonode, nhdplusreg fields
        Raster File (tif) - user must specify if the data is continuous or classed 
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

"""
Method: readData()
Purpose: Read all necessary input data and return it. Read in the NHD+ catchments, raster, and
         routing database
Params: NHDpath - path to NHD+ catchments (shapefile)
        sassDB - path to routing database (SASS database or CSV)
        sassType - flag for routing database format - 0 for SASS db and 1 for CSV
        rasPath - path to raster
        rasType - flag for raster type - 0 for continuous, 1 for classed
        noData - raster no data value
        VPUzone - VPU zone (from NHD) that the catchments exist in. Used to shrink 
                  routing database from CONUS scale.
Returns: NHD - geodataframe of NHD+ catchments
         shapes - dictionary of fiona geometries for NHD+ catchments
         db - pandas dataframe of routing database
         rasVals - list of unique raster values for classed data (empty list for continuous)
"""
def readData(NHDpath, sassDB, sassType, rasPath, rasType, noData, VPUzone):
    try:
        print("Reading Data...")
        # Read catchments as geodataframe -- used to join with final results
        NHD = gpd.read_file(NHDpath) #read shapefile as geopandas df
        NHD = NHD[['FEATUREID', 'geometry']]
        # Read routing database 
        if sassType == 0:
            with SAS7BDAT (sassDB, skip_header=False) as r:
                db = r.to_data_frame()
                db = db[['FL_ComID', 'cfromnode', 'ctonode', 'nhdplusreg']]
                db = db[db.nhdplusreg == VPUzone]
                db.to_csv("/data/gis2/sm/NHD_Test/Version3_MWDB/Output/SassDB_02.csv", index=False)
                db = db.fillna(-100)
                db["cfromnode"] = db.cfromnode.astype(int)
                db["ctonode"] = db.ctonode.astype(int)
        elif sassType == 1:
            db = pd.read_csv(sassDB)
            db = db.dropna(axis='index')
            db['FL_ComID'] = db['FL_ComID'].apply(int)
        else:
            print("Invalid sassType\n\t0 for .sas7bdat file\n\t1 for .csv file")
        # Read catchments as fiona geometries and store in dictionary
        with fiona.open(NHDpath, "r") as geoms:
            shapes = {feature['properties']['FEATUREID']:feature['geometry'] for feature in geoms}
        # If the raster is classed -- get all unique class values and store in list
        if rasType == 1: #discrete data -- get unique values
            ras = gdal.Open(rasPath)
            rasVals = list(np.unique(np.array(ras.GetRasterBand(1).ReadAsArray())))
            if noData in rasVals:
                rasVals.remove(noData)
        else:
            rasVals = []
        #return data
        return NHD, shapes, db, rasVals
    except:
        print("Error Reading Input Data -- Check Paths\n")
        return pd.DataFrame(), {}, pd.DataFrame(), []

"""
Method: aggregate()
Purpose: calculate raster statistics for each catchment, find the upstream network for each catchment and
         accumulate the raster statistics.
Params: shapes - dict of fiona geometries for the catchments
        PLUS - dataframe containing the routing database
        stat - statistic type to be calculated for the catchments (only needed for continuous data)
        rasPath - path to the raster data
        convFacct - the conversion factor specified by the user
        noData - the raster noData value
        rasType - 0 or 1 specifying if raster is continuous or classed
        rasVals - list containing unique class values for a classed raster (empty list for continuous)
Returns: sumDict - dictionary where the catchment ID is the key and the data is a list of catchment statistics and
                   accumulated statistics
"""
def aggregate(shapes, PLUS, stat, rasPath, convFact, noData, rasType, rasVals):
    print("Calculating...")
    sumDict = {}
    #Calculate statistic for each catchment
    with rio.open(rasPath) as src:
        if rasType == 0: #continuous raster
            for s in shapes:
                val, ar = getSum(src, shapes[s], stat, convFact, noData)
                sumDict[s] = [float(val), float(ar)]
        elif rasType == 1: #classed raster
            for s in shapes:
                vals = getSumDis(src, shapes[s], noData, rasVals)
                sumDict[s] = vals
                
    #Create dict of direct upstream connections -- faster than referencing table directly each time
    allUp = {}
    allLines = list(PLUS['FL_ComID'])
    for s in allLines:
        fromN = list(PLUS[PLUS['FL_ComID'] == s]['cfromnode'])[0]
        if -100 == fromN: #Nan values converted to -100; ignore these
            allUp[s] = []
        else:
            allUp[s] = list(PLUS[PLUS['ctonode'] == fromN]['FL_ComID'])

    #Accumulate statistic for each catchment
    print("Accumulating...")
    allCatch = list(shapes.keys())
    if rasType == 0: #if continuous data -- update list of "unique values" to be the statistic and area
        rasVals = ['Stat', 'Area']
    numClasses = len(rasVals)
    for s in shapes:
        upstream = allUpstream(allUp, s, PLUS) #all upstream connections
        upCatch = list((set(upstream)&set(allCatch))) #all upstream catchments
        sumDict[s] = getSumList(sumDict, upCatch, numClasses, convFact, s) #list of catchment values (original scale) and accumulated values (converted scale)
    return sumDict #return all data in dict

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
                FL_ComID in the routing database
        s - the current catchment ID
        PLUS - dataframe of the routing database
Returns: toSum - list of complete upstream connections (catchments and flowlines)
"""
def allUpstream(allUp, s, PLUS):
    toSum = [int(s)] #was allUp[s] + [s] changed since moving allUp to iterCheck
    new = [int(s)]
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
Purpose: find direct upstream connections for all FL_ComIDs in a list and build a list
         of the unique IDs found. Called by allUpstream
Params: coms - list of FL_ComIDs to find direct upstream connections of
        allUp - dictionary containing the direct upstream connections for each unique
                FL_ComID in the routing database
        PLUS - dataframe of the routing database
Returns: new - list of unique upstream catchment IDs
         notFound - list of unique IDs that did not exist in allUp
"""
def iterCheck(coms, allUp, PLUS):
    new = []
    notFound = []
    for c in coms:
        if c not in allUp:
            fromN = list(PLUS[PLUS['FL_ComID'] == c]['cfromnode'])[0]
            if -100 == fromN: #Nan values converted to -100; ignore these
                allUp[c] = []
            else:
                allUp[c] = list(PLUS[PLUS['ctonode'] == fromN]['FL_ComID'])
                allUp[c] = [int(x) for x in allUp[c]]
        if c in allUp:
            if len(allUp[c]) > 0:
                new = new + allUp[c]
        else:
            print(c, ": not in allUp")
            notFound.append(c)
    return list(set(new)), notFound

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
        print("Could Not Mask Raster -- Check Projections of Files")
        return []

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
def getSum(src, geom, statType, convFact, noData):
    try:
        ary, t = rio.mask.mask(src, [geom], crop=True, all_touched=False) #mask by current catchment
        area = ((ary != noData).sum())* convFact # number of pixels that are not no data
        #set no data value to nan for float rasters and 0 for integer rasters
        if type(ary[0][0][0]) == np.float32 or type(ary[0][0][0]) == np.float64:
            ary[ary == noData] = np.nan #set nodata to nan 
        else:
            ary[ary == noData] = 0
        if statType == "MAX":
            catchSum = float(np.amax(ary))  #sum ras in current catchment -- CHANGE BACK TO SUM
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
        print("Could Not Mask Raster -- Check Projections of Files")
    return 0

"""
Method: writeTable()
Purpose: Produce the catchment statistics and accumulated statistics as a CSV file
Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
        path - directory to write the CSV to
        name - desired filename (must end with .csv)
        rasVals - list of unique values for classed data (used as column names)
Returns: None
"""
def writeTable(sumDict, path, name, rasVals):
    print("Writing CSV...")
    if name[-4:] == ".csv":
        if len(rasVals) == 0:
            rasVals = ['Stat', 'Area']
        cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
        for idx, v in enumerate(rasVals):
            cols[idx] = str(rasVals[idx])
            cols[idx+len(rasVals)] = "Acc_"+str(rasVals[idx])
        df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
        df['FEATUREID'] = df.index
        df = df.sort_values(by=['FEATUREID'])
        finCols = ['FEATUREID'] + cols
        df = df[finCols]
        fil = os.path.join(path, name)
        df.to_csv(fil, index=False)
    else:
        print("Output File Must be of Type .csv")

"""
Method: writeShapefile()
Purpose: Produce the catchment statistics and accumulated statistics as a shapefile. Joins
         the statistics with the NHD+ catchments.
Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
        NHD - geopandas dataframe for the catchments
        path - directory to write the CSV to
        name - desired filename (must end with .csv)
        rasVals - list of unique values for classed data (used as column names)
Returns: None
"""
def writeShapefile(sumDict, NHD, path, name, rasVals):
    print("Writing Shapefile...")
    if name[-4:] == ".shp":
        f = os.path.join(path, name)
        if len(rasVals) == 0: #continous data
            rasVals = ['Stat', 'Area']
        cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
        for idx, v in enumerate(rasVals):
            cols[idx] = str(rasVals[idx])
            cols[idx+len(rasVals)] = "Acc_"+str(rasVals[idx])
        df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
        df['FEATUREID'] = df.index
        NHD = NHD.merge(df, on=['FEATUREID'], how='left')
        NHD.to_file(f)
    else:
        print("Output File Must be of Type .shp")

    
"""
Global variables needing user input
rasType: 0 (continuous raster) or 1 (classified raster)
statType: "MIN", "MAX", "MEDIAN", "MEAN", or "SUM" (choose one statistic for continuous raster only)
outputType: [1, 1] (first index is output type CSV and second is shapefile; 0 is do NOT write and 1 is write)
cellSize: cell size of the raster (will automate this later)
noData: no data value of the raster (will automate this later)
sassType: 0 (routing database is in sass db .sas7bdat) or 1 (routing database is in CSV file)
convFact: conversion factor to apply to accumulated values; recommend output be in sq km
VPUzone: NHD zone that the catchments are in (02 is MidAtlantic)
"""
rasType = 0 #0 for continuous; 1 for classed
#CAN DO THESE : min (amin), max (amax), median (median), mean (mean)
statType = "SUM" # -- only useful for continuous data
#NEED TO KNOW HOW TO PRODUCE DATA -- CSV, SHAPEFILE, or BOTH
outputType = [1, 1]
cellSize = 10
noData = 65535 #-128 
sassType = 1 #1 is csv, 0 is .sas7bdat
convFact = (cellSize * cellSize) / 1000000 #conversion from sq meter to sq km #0.00404686 #acres to sq km 
VPUzone = "02" #Mid atlantic (ches bay and delaware)

"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    #Directory Paths
    ABS_DIR = os.path.abspath("/data/gis2/sm/NHD_Test/Version3_MWDB") #Main Directory 
    IN_DIR = os.path.join(ABS_DIR, "Input") #Directory containing the input data (catchments, routing db and raster)

    #File Paths
    NHD_SHP = os.path.join(IN_DIR, "Potomac_Catchmets_P.shp") 
    sassPath = os.path.join(ABS_DIR, "Output", "SassDB_02.csv") #nhdplusv2_us.sas7bdat
    RAS = os.path.join(IN_DIR, "IMP_P6LU.tif") #input raster

    #Read input data
    startTime = float(time.time())
    NHD, shapes, PLUS, rasVals = readData(NHD_SHP, sassPath,  sassType, RAS, rasType, noData, VPUzone)
    lastTime = float(time.time())
    print("\nTime: ", str((lastTime- startTime)/60), " min\n")

    #Calculate statistics and accumulate
    data = aggregate(shapes, PLUS, statType, RAS, convFact, noData, rasType, rasVals)
    print("\nTime: ", str((float(time.time()) - lastTime)/60), " min\n")

    #write out data
    lastTime = float(time.time())
    if outputType[0] == 1:
        writeTable(data, os.path.join(ABS_DIR, "Output"), "POT_IMP_81320.csv", rasVals)
    if outputType[1] == 1:
        writeShapefile(data, NHD, os.path.join(ABS_DIR, "Output"), "POT_IMP_81320.shp", rasVals)
    print("\nTotal Time: ", str((lastTime - startTime)/60), " min\n")
   
    
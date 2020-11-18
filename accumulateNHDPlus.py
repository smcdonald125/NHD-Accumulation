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

####################################################################################################
####################       READ INPUT DATA       ###################################################
####################################################################################################

"""
Method: readCatchmentGeoms()
Purpose: Read NHD+ catchment geometries into a dictionary
Params: NHDpath - path to NHD+ catchments (shapefile)
Returns: shapes - dictionary of fiona geometries for NHD+ catchments
"""
def readCatchmentGeoms(NHDpath, catchID):
    print("Reading catchment geometries...")
    if os.path.isfile(NHDpath):
        with fiona.open(NHDpath, "r") as geoms:
            shapes = {feature['properties'][catchID]:feature['geometry'] for feature in geoms}
        return shapes
    else:
        print("Invalid File Path: ", NHDpath)
        print("Could not read catchments shapefile")
        return {}

"""
Method: readRouting()
Purpose: Read routing database into a dataframe and return it.
Params: sassDB - path to routing database (SASS database or CSV)
        VPUzone - VPU zone (from NHD) that the catchments exist in. Used to shrink 
                  routing database from CONUS scale.
Returns: db - pandas dataframe of routing database
"""
def readRouting(sassDB, VPUzone, NHD_lines):
    print("Reading routing database...")
    if os.path.isfile(sassDB):
        if sassDB[-8:] == 'sas7bdat':
            db = pd.read_sas(sassDB)
            db = db[['FL_ComID', 'cfromnode', 'ctonode', 'nhdplusreg', 'DivFrac', 'ModDivFrac', 'div2ground', 'divgroup', 'Hydroseq']]
            db['nhdplusreg'] = db['nhdplusreg'].str.decode("utf-8")
            db = db[db.nhdplusreg == VPUzone]
            db = db.fillna(-100)
            db["cfromnode"] = db.cfromnode.astype(int)
            db["ctonode"] = db.ctonode.astype(int)
            lines = gpd.read_file(NHD_lines)
            lines = lines[['COMID', 'FTYPE']]
            coast = list(lines[lines['FTYPE'] == 'Coastline']['COMID'])
            downNodes = list(db[db['FL_ComID'].isin(coast)]['ctonode']) #list of nodes coming from coastline
            db.loc[db['cfromnode'].isin(downNodes), 'cfromnode'] = -100 #can't remove recrods -- update the from to be null
            db.to_csv("/data/gis2/sm/NHD_Test/Version3_MWDB/Output/SassDB_"+str(VPUzone)+".csv", index=False)
        elif sassDB[-3:] == 'csv':
            db = pd.read_csv(sassDB)
            db = db.dropna(axis='index')
            db['FL_ComID'] = db['FL_ComID'].apply(float)
            db["cfromnode"] = db.cfromnode.astype(int)
            db["ctonode"] = db.ctonode.astype(int)
        else:
            print("Invalid sassType\n\t0 for .sas7bdat file\n\t1 for .csv file")
            return pd.DataFrame()
        return db
    else:
        print("Invalid File Path: ", sassDB)
        return pd.DataFrame()

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
        ras = gdal.Open(rasPath)
        gt = ras.GetGeoTransform() 
        cellSize = gt[1] 
        if rasType == 1:
            rasVals = list(np.unique(np.array(ras.GetRasterBand(1).ReadAsArray())))
            return cellSize, rasVals
        else:
            return cellSize, []
    else:
        print("Invalid file path: ", rasPath)
        return 0, []

"""
Method: readTabular()
Purpose: create data frame with catchment statistics from an existing table
Params: tabPath - path to the raster data
        idName - string of column name containing the unique COMIDs
        cols - list of columns names to accumulate
Returns: sumDict - dictionary where the catchment ID is the key and the data is a list of catchment statistics and
                   accumulated statistics
"""
def readTabular(tabPath, idName, cols):
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
####################       READ / CALCULATE CATCHMENT STATISTICS       #############################
####################################################################################################

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
def aggregateRas(shapes, stat, rasPath, cellSize, rasType, rasVals):
    print("Calculating Catchment Statistics...")
    sumDict = {}
    #Calculate statistic for each catchment
    with rio.open(rasPath) as src:
        noData = src.nodatavals
        if noData in rasVals: #ignore NoData
            rasVals.remove(noData)
        if rasType == 0: #continuous raster
            for s in shapes:
                val, ar = getSum(src, shapes[s], stat, cellSize, noData)
                sumDict[s] = [float(val), float(ar)]
        elif rasType == 1: #classed raster
            for s in shapes:
                vals = getSumDis(src, shapes[s], noData, rasVals)
                sumDict[s] = vals
    return rasVals, sumDict


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
        cellSize - raster cell size to multiply by total pixel count
        noData - raster noData value 
Returns: catchSum - catchment statistic
         area - catchment area
"""
def getSum(src, geom, statType, cellSize, noData):
    try:
        ary, t = rio.mask.mask(src, [geom], crop=True, all_touched=False) #mask by current catchment
        area = ((ary != noData).sum())*(cellSize*cellSize) # number of pixels that are not no data - removed * convFact 
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
        
####################################################################################################
####################       ACCUMULATE DATA       ###################################################
####################################################################################################
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

####################################################################################################
####################       GENERATE UPSTREAM NETWORK      ##########################################
####################################################################################################

"""
Method: generateUpstream()
Purpose: To compute the complete upstream network for the specified catchment
Params:  PLUS - dataframe of the routing database
Returns: toSum - list of complete upstream connections (catchments and flowlines)
"""
def generateUpstream(PLUS):
    print("Finding direct upstream links...")
    #Create dict of direct upstream connections -- faster than referencing table directly each time
    allUp = {}
    allLines = list(set(list(PLUS['FL_ComID']))) #get all unique IDs in relation table
    for s in allLines:
        fromN = list(PLUS[PLUS['FL_ComID'] == s]['cfromnode'])[0]
        if -100 == fromN: #Nan values converted to -100; ignore these
            allUp[s] = []
        else:
            allUp[s] = list(PLUS[PLUS['ctonode'] == fromN]['FL_ComID'])
    return allUp

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
    toSum = [s] #was allUp[s] + [s] changed since moving allUp to iterCheck
    new = [s]
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
                #allUp[c] = [int(x) for x in allUp[c]]
        if c in allUp:
            if len(allUp[c]) > 0:
                new = new + allUp[c]
        else:
            print(c, ": not in allUp")
            notFound.append(c)
    return list(set(new)), notFound


####################################################################################################
#########################       WRITE DATA       ###################################################
####################################################################################################

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
        print("Output File Must be of Type .csv")

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
            NHD = NHD.merge(df, on=['COMID'], how='left')
            NHD.to_file(f)
        else:
            print("Output File Must be of Type .shp")
    else:
        print("Invalid File Path: ", NHDPath)
        print("No catchment shapefile to join with data")



    
"""
Global variables needing user input
rasType: 0 (continuous raster) or 1 (classified raster) or 2 (tabular)
statType: "MIN", "MAX", "MEDIAN", "MEAN", or "SUM" (choose one statistic for continuous raster only)
outputType: [1, 1] (first index is output type CSV and second is shapefile; 0 is do NOT write and 1 is write)
convFact: conversion factor to apply to accumulated values; recommend output be in sq km
VPUzone: NHD zone that the catchments are in (02 is MidAtlantic)
catchID: Name of column with unique catchment ID in the catchments shapefile
tabID: if using tabular data containing catchment statistics, this is the name of the column with the catchment IDs
tabCols: if using tabular data containing catchment statistics, this is a list of columns to accumulate
"""
rasType = 2 #0 for continuous; 1 for classed; 2 fo tabular
statType = "SUM" # -- only useful for continuous data
outputType = [0, 1] #what format to produce data
convFact = 1 #acres to sq km  #1 / 1000000 #conversion from sq meter to sq km 
VPUzone = "02" #Mid atlantic (ches bay and delaware)
catchID = 'COMID'
tabID = 'COMID' #only need for rasType 2
tabCols = ['CRP', 'FORE', 'INR', 'IR', 'MO', 'PAS', 'TCI', 'TCT', 'TG', 'WAT', 'WLF', 'WLO', 'WLT'] #only need for rasType 2


"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    invalidData = False #Flag to clean exit if any function fails
    #Directory Paths
    ABS_DIR = os.path.abspath("/data/gis2/sm/NHD_Test/Version3_MWDB") #Main Directory 
    OUT_DIR = r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/AccumulateNHDPlus/P6LU2013"
    IN_DIR = os.path.join(ABS_DIR, "Input") #Directory containing the input data (catchments, routing db and raster)

    #File Paths
    NHD_SHP = os.path.join(IN_DIR, "NHDv2_CHWAv2.shp") 
 
    sassPath = os.path.join(IN_DIR, "nhdplusv2_us.sas7bdat")  #"SassDB_02.csv") 
    NHD_lines = os.path.join(IN_DIR, "NHDv2_CHWALines_P.shp") #NEED TO REMOVE COASTLINE CONNECTIONS -- only need to do this if passing sass db -- csv produed will have this done
    inputData = os.path.join(IN_DIR, "IMP_P6LU.tif") # "NHDv2_CHWAv2.shp")#input raster or tabular catchment data
    
    #Read input data
    startTime = float(time.time())
    sumDict = {}
    upstream = {} #build once - code loop to run multiple datasets
    PLUS = readRouting(sassPath, VPUzone, NHD_lines)
    if len(PLUS) == 0:
        invalidData = True
    #Generate upstream network
    if not invalidData:
        allUp = generateUpstream(PLUS)
        if len(allUp) == 0:
            invalidData = True
    if not invalidData:     
        if rasType == 0 or rasType == 1: #Raster data - continuous or categorical
            cellSize, rasVals = readRas(inputData, rasType)
            if cellSize == 0:
                invalidData = True
            else:
                shapes = readCatchmentGeoms(NHD_SHP, catchID)
                if len(shapes) == 0:
                    invalidData = True
                else:
                    rasVals, sumDict = aggregateRas(shapes, statType, inputData, cellSize, rasType, rasVals)
                    if len(sumDict) == 0:
                        invalidData = True
                    else:
                        sumDict = accumulate(rasType, rasVals, allUp, sumDict, convFact)
        elif rasType == 2: #tabular data
            sumDict = readTabular(inputData, tabID, tabCols)
            if len(sumDict) == 0:
                invalidData = True
            else:
                sumDict, upstream = accumulate(rasType, tabCols, allUp, upstream, sumDict, convFact) 
                if len(sumDict) == 0:
                    invalidData = True
                        

    lastTime = float(time.time())
    print("\nTime: ", str((lastTime- startTime)/60), " min\n")

    #write out data
    if rasType == 2:
        rasVals = tabCols.copy()
    if not invalidData:
        if outputType[0] == 1:
            writeTable(sumDict, OUT_DIR, "CBW_IMP_090420.csv", rasVals)
        if outputType[1] == 1:
            writeShapefile(sumDict, NHD_SHP, OUT_DIR, "CBW_P6LU2013_103020.shp", rasVals)
    print("\nTotal Time: ", str((lastTime - startTime)/60), " min\n")
   
    

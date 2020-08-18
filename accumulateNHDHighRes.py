"""
Purpose: Calculate raster statistics at the NHD+ catchment level and accumulate the data. The input raster can be continuous data,
        where the following statistics can be calculated: MAX, MIN, MEAN, MEDIAN, SUM. The raster can also be classed data, where the
        total area per class is calculated. The catchment statistics and accumulated data can be written as a table to a CSV and
        can be written as a shapefile, joined to the NHD+ catchments. The upstream networks are built using a routing database, produced
        by Mike Wieczorek. The routing database can be reduced to the VPU zone level, and can be read in as its original form (SASS database)
        or a CSV file. 
Data:   NHD High Res catchments as shapefile - Shapefile to act as zones
        PlusFlow CSV - NHD PlusFlow table to build the upstream networks
                        OR
        ecoSheds high res catchments as shapefile - Shapefile to act as zones; Can be used to build upstream network
        ecodSheds high res flowline dbf as CSV - can be used to build upstream network (optional if using ecoSheds catchments)

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
from osgeo import gdal
import fiona
import math
import os
import time

"""
Method: readData24k()
Purpose: Read all necessary input data and return it. Read in the high res catchments, raster, and
         table with connectivty (nhd high res - PlusFlow, ecoSheds - catchments dbf or flowlines dbf as CSV)
Params: Catchpath - path to high res catchments (shapefile)
        plusTabPath - path to table with connectivity data (CSV)
        interVPUPath - Table connecting HUC4 zones for NHD High Res data (Created manually - included in NHD PlusFlow?)
        VPUzone - HUC4 zone ID; Reduce PlusFlow to HUC4 zone (for nhd high res only)
        rasPath - path to raster
        rasType - flag for raster type - 0 for continuous, 1 for classed
        noData - raster no data value
        catchName - name of the unique catchment ID in the catchments shapefile
        isECO - boolean flag to denote if using NHD high res data or ecoSHEDs data
Returns: catch - geodataframe of NHD+ catchments
         PLUS - pandas dataframe that contains the connectivity data
         shapes - dictionary of fiona geometries for NHD+ catchments
         rasVals - list of unique raster values for classed data (empty list for continuous)
"""        
def readData24k(Catchpath, plusTabPath, interVPUPath, VPUzone, rasPath, rasType, noData, catchName, isECO):
    try:
        print("Reading data...")
        #read catchments as geodataframe
        catch = gpd.read_file(Catchpath) 
        catch = catch[[catchName, 'geometry']] 
        PLUS = pd.read_csv(plusTabPath) 
        if not isECO:
            catch = catch.dropna(axis='index')
            catch['TMPNAME'] = catch[catchName].apply(float)
            catch = catch[['TMPNAME', 'geometry']]
            catch = catch.rename(columns={'TMPNAME':catchName})
            PLUS = PLUS[['FromNHDPID', 'ToNHDPID', 'FromVPUID', 'ToVPUID']]
            #interVPU = pd.read_csv(interVPUPath)
            #PLUS = PLUS.append(interVPU, sort=False)
            if VPUzone > 0:
                print("PLUS reduced to zone: ", VPUzone)
                PLUS = PLUS[PLUS['ToVPUID'] == float(VPUzone)]
                PLUS = PLUS[(PLUS['FromVPUID'] == float(VPUzone)) | (PLUS['FromNHDPID'] == float(0))]
        else:
            PLUS = PLUS[['FEATUREID', 'NextDownID']]
        
        #shapely geoms for catch
        shapes = {}
        with fiona.open(Catchpath, "r") as geoms:
            shapes = {feature['properties'][catchName]:feature['geometry'] for feature in geoms} #dict of geomtries with FeatureID as key
       
        # If the raster is classed -- get all unique class values and store in list
        if rasType == 1: #discrete data -- get unique values
            ras = gdal.Open(rasPath)
            rasVals = list(np.unique(np.array(ras.GetRasterBand(1).ReadAsArray())))
            if noData in rasVals:
                rasVals.remove(noData)
        else:
            rasVals = []

        return catch, PLUS, shapes, rasVals
    except:
        print("Error Reading Input Data -- Check Paths\n")
        return pd.DataFrame(), pd.DataFrame(), {}, []
       
"""
Method: calcSum()
Purpose: calculate raster statistics for each catchment, find the upstream network for each catchment and
         accumulate the raster statistics.
Params: shapes - dictionary of fiona geometries for the catchments
        PLUS1 - pandas dataframe with connectivity info
        RasPath - path to raster
        convFact - conversion factor to be applied to accumulated data
        noData - no data value of the raster
        toFromNames - list of length 2 containing the names of the to/from columns; [to, from]
        isECO - boolean flag denoting if using NHD high res or ecoSHEDs
        rasVals - list of unique raster values of classified raster (empty list for continuous data)
"""
def calcSum(shapes, PLUS1, RasPath, convFact, noData, toFromNames, isECO, rasVals):
    #Get list of "headwaters" -- starts of flow paths
    hw = list(PLUS1[PLUS1[toFromNames[1]] == 0][toFromNames[0]]) #[1] = from [0] = to
    hw = [float(x) for x in hw]
    NHD_HW = list(shapes.keys())
    NHD_HW = [float(x) for x in NHD_HW if x != None]
    if 0 in NHD_HW:
        NHD_HW.remove(0)
    no_hw = list(set(NHD_HW) - set(hw)) 
    #Dictionary to store sum, aggregated sum with FeatureID as key -- easy to convert to DF and merge
    sumDict = {}
    allUp = {}
    print("Calculating Statistic...")
    global statType
    global rasType
    #Mask the raster by each polygon in a loop and calculate sum
    with rio.open(RasPath) as src:
        for cID in NHD_HW:
            if rasType == 0: #continuous raster
                if not isECO:   
                    val, a = getSum(src, [shapes[str(int(cID))]], statType, convFact, noData) 
                else:
                    val, a = getSum(src, [shapes[cID]], statType, convFact, noData)
                sumDict[float(cID)] = [val, a]
            elif rasType == 1: #calssified raster
                if not isECO:
                    vals = getSumDis(src, shapes[str(cID)], noData, rasVals)
                else:
                    vals = getSumDis(src, shapes[cID], noData, rasVals)
                sumDict[float(cID)] = vals
            #find direct upstream connections    
            x = list(PLUS1[PLUS1[toFromNames[0]] == float(cID)][toFromNames[1]])
            x = [float(i) for i in x]
            x = list(set(x))
            if 0 in x:
                x.remove(0)
            allUp[float(cID)] = x

    #add line IDs in PLUS that are not catchments to allUp -- this may happen for 24k data
    allPIDs = list(set(list(PLUS[toFromNames[0]]) + list(PLUS[toFromNames[1]]))) # all IDs in the PlusFlow (NHD) Catchments/flowlines dbf (ecosheds)
    allPIDs = [float(x) for x in allPIDs]
    if 0 in allPIDs:
        allPIDs.remove(0)
    toAdd = list(set(allPIDs) - set(list(allUp))) # IDs that are connected but aren't yet in the upstream dictionary
    for ta in toAdd: #add missing data to the dictionary
        x = list(PLUS1[PLUS1[toFromNames[0]] == ta][toFromNames[1]])
        x = [float(i) for i in x]
        x = list(set(x))
        if 0 in x:
            x.remove(0)    
        allUp[float(ta)] = x

    print("Accumulating Statistic...")
    allCatch = list(shapes.keys())
    if rasType == 0: #if continuous data -- update list of "unique values" to be the statistic and area
        rasVals = ['Stat', 'Area']
    numClasses = len(rasVals)
    for cID in no_hw:
        upstream = allUpstream(allUp, cID, PLUS, isECO) #all upstream connections
        upCatch = list((set(upstream)&set(allCatch))) #all upstream catchments
        sumDict[float(cID)] = getSumList(sumDict, upCatch, numClasses, convFact, cID) #list of catchment values (original scale) and accumulated values (converted scale)        
    
    return sumDict


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
        isECO - boolean flag denoting if the high res data is ecoSHEDs or NHD high res
        
Returns: toSum - list of complete upstream connections (catchments and flowlines)
"""
def allUpstream(allUp, s, PLUS, isECO):
    toSum = [float(s)] #was allUp[s] + [s] changed since moving allUp to iterCheck
    new = [float(s)]
    while len(new) > 0:
        new, x = iterCheck(new, allUp, PLUS, isECO) #iterative version of recursive check()
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
        isECO - boolean flag denoting if the high res data is ecoSHEDs or NHD high res
Returns: new - list of unique upstream catchment IDs
         notFound - list of unique IDs that did not exist in allUp
"""
def iterCheck(coms, allUp, PLUS, isECO):
    new = []
    notFound = []
    if isECO:
        fromName = 'FEATUREID'
        toName = 'NextDownID'
    else:
        fromName = 'FromNHDPID'
        toName = 'ToNHDPID'
    for c in coms:
        if c not in allUp:
            x = list(PLUS[PLUS[toName] == c][fromName])
            x = [float(y) for y in x]
            allUp[c] = x
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
def getSum(src, geom, stat, convFact, noData):
    try:
        ary, t = rio.mask.mask(src, geom, crop=True, all_touched=False) #mask by current catchment
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
        print("Could Not Mask Raster -- Check Projections of Files")
    return 0
    

"""
Method: writeTable()
Purpose: Produce the catchment statistics and accumulated statistics as a CSV file
Params: sumDict - dictionary with all statistics and accumulated statistics for the catchments
        path - directory to write the CSV to
        name - desired filename (must end with .csv)
        rasVals - list of unique values for classed data (used as column names)
        catchIDName - string denoting the catchment ID field name of the high res catchments
Returns: None
"""
def writeTable(sumDict, path, name, rasVals, catchIDName):
    print("Writing CSV...")
    if name[-4:] == ".csv":
        if len(rasVals) == 0:
            rasVals = ['Stat', 'Area']
        cols = [0] * (len(rasVals) * 2) #list to store catchment stat and acc stats
        for idx, v in enumerate(rasVals):
            cols[idx] = str(rasVals[idx])
            cols[idx+len(rasVals)] = "Acc_"+str(rasVals[idx])
        df = pd.DataFrame.from_dict(sumDict, orient='index', columns=cols)
        df[catchIDName] = df.index
        df = df.sort_values(by=[catchIDName])
        finCols = [catchIDName] + cols
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
        catchIDName - string denoting the catchment ID field name of the high res catchments
Returns: None
"""
def writeShapefile(sumDict, NHD, path, name, rasVals, catchIDName):
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
        df[catchIDName] = df.index
        NHD = NHD.merge(df, on=[catchIDName], how='left')
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
convFact: conversion factor to apply to accumulated values; recommend output be in sq km
catchIDName: name of the columns denoting the unique catchment IDs
toFromPlusNames: list of length to denoting the column names for the to/from columns;
                 for NHD high res data ['ToNHDPID', 'FromNHDPID']
                 for ecoSHEDs ['NextDownID', 'FEATUREID'] 
VPUzone: NHD HUC4 zone that the catchments are in (207 is Potomac);
         0 if larger than one zone
isECO: Boolean flag to denote which version of catchments/connections are being used.
       True - using ecoSHEDs data
       False - using NHD High Res catchments and PlusFlow
"""
rasType = 0
statType = "SUM"
outputType = [1, 1]
cellSize = 10
noData = 65535
convFact = (cellSize * cellSize) / 1000000.0 
catchIDName = 'COMID'
toFromPlusNames = ['ToNHDPID', 'FromNHDPID'] 
VPUzone = 207  
isECO = False

"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    #Directory Paths
    ABS_DIR = os.path.abspath("/data/gis2/sm/NHD_Test/Version2_SM/24k") 
    IN_DIR = os.path.join(ABS_DIR, "Input")

    #Data Paths
    NHD_SHP = os.path.join(IN_DIR, "Catchments_207_Proj.shp") #Path to Catchments shapefile
    PLUS_TAB = os.path.join(IN_DIR, "PlusFlow_207.csv") #Path to connecitivity data (ecoSheds export catchment or flowlines dbf to csv)
    RAS = os.path.join(IN_DIR, "IMP_P6LU.tif") #path to raster
    interVPUP = os.path.join(IN_DIR, "interVPU_HUC4.csv") #path to interVPU table; can be ignored if within one HUC4; Not needed for ecoSHEDS
    
    #Read data
    startTime = float(time.time())        
    catch, PLUS, shapes, rasVals = readData24k(NHD_SHP, PLUS_TAB, interVPUP, VPUzone, RAS, rasType, noData, catchIDName, isECO)
    lastTime = float(time.time())
    print("\nTime: ", str((lastTime- startTime)/60), " min\n")

    #Calculate catchments statistis and accumulate
    sumDict = calcSum(shapes, PLUS, RAS, convFact, noData, toFromPlusNames, isECO, rasVals) # -- sumDis is only for discrete data
    print("\nTime: ", str((lastTime- time.time())/60), " min\n")
    lastTime = float(time.time())
    
    #Write data
    if outputType[0] == 1:
        writeTable(sumDict, os.path.join(ABS_DIR, "Output"), "24k_POT_081120.csv", rasVals, catchIDName)
    if outputType[1] == 1:
        writeShapefile(sumDict, catch, os.path.join(ABS_DIR, "Output"), "24k_POT_081120.shp", rasVals, catchIDName)
    print("\nTime: ", str((lastTime- time.time())/60), " min\n")
    lastTime = float(time.time())
    print("Total Time: ", str((lastTime - startTime)/60), " min")
  

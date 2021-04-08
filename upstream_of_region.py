"""
Purpose: Find the all catchments draining to a specified region. The upstream networks are built using a routing database, produced
        by Mike Wieczorek. The routing database can be reduced to the VPU zone level, and can be read in as its original form (SASS database)
        or a CSV file. 
Data:   NHD+ catchments as shapefile - Shapefile to act as zones - must have COMID
        Routing Database as SASS database or CSV- Must have FL_ComID, cfromnode, ctonode, nhdplusreg fields
        Catchments in region of interest shapefile - drainage zone
Authors:  Sarah McDonald, Geographer, USGS
Contact: smcdonald@chesapeakebay.net
"""

#Import Libs
import pandas as pd
import geopandas as gpd
import os
from sas7bdat import SAS7BDAT
import multiprocessing as mp

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
    print("Reading scatchment geometries...")
    if os.path.isfile(NHDpath):
        with fiona.open(NHDpath, "r") as geoms:
            shapes = {feature['properties'][catchID]:feature['geometry'] for feature in geoms}
        return shapes
    else:
        print("Invalid File Path: ", NHDpath)
        print("Could not read catchments shapefile")
        return {}

def readRouting(sassDB, VPUzone, NHD_lines):
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
            db.to_csv(r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/ROUTING_DATA/SassDb_02_040721.csv", index=False)
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
        
####################################################################################################
####################       GENERATE UPSTREAM NETWORK      ##########################################
####################################################################################################

def getUpstreamMP(allUp, regionCOMIDs, PLUS):
    """
    Method: getUpstreamMP()
    Purpose: Chunk the comids in the region and find their upstream catchments on multiple cores.
    Params: allUp - dictionary of direct upstream links
            regionCOMIDs - list of COMIDs in the region
            PLUS - df of routing db
    Returns: fin - list of COMIDs that drain to catchments in region
    """
    #pull cores
    numCores = mp.cpu_count() - 2

    #chunk up data
    batch_size = int((len(regionCOMIDs) / numCores) + 1)
    chunk_iterator = []
    for i in range(numCores):
        mn, mx = i * batch_size, (i + 1) * batch_size
        args = regionCOMIDs[mn:mx], allUp, PLUS
        chunk_iterator.append(args)
    
    #call intersect function
    pool = mp.Pool(processes=numCores)
    results = pool.map(allUpstream, chunk_iterator)
    pool.close()

    #concat results
    fin = []
    for r in results:
        fin += r

    return fin

def allUpstream(args):
    """
    Method: allUpstream()
    Purpose: To find all upstream catchments in the list
    Params: args - tuple of:
                catchmnts - list of catchment IDs
                allUp - dictionary containing the direct upstream connections for each unique
                    FL_ComID in the routing database
                PLUS - dataframe of the routing database
    Returns: toSum - list of complete upstream connections (catchments and flowlines)
    """
    catchmnts, allUp, PLUS = args
    allCatch = []
    for s in catchmnts:
        toSum = [int(s)] 
        new = [int(s)]
        while len(new) > 0:
            new, x = iterCheck(new, allUp, PLUS) 
            if len(x) > 0:
                toSum = list(set(toSum) - set(x))
                new = list(set(new) - set(x))
            new = list(set(new) - set(toSum)) #remove any items from new that are already stored in the network
            toSum = toSum + new
        allCatch += list(set(toSum))
    return list(set(allCatch))

def generateUpstream(PLUS):
    """
    Method: generateUpstream()
    Purpose: To compute the complete upstream network for the specified catchment
    Params:  PLUS - dataframe of the routing database
    Returns: toSum - list of complete upstream connections (catchments and flowlines)
    """
    print("Find direct upstream links...")
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

def iterCheck(coms, allUp, PLUS):
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

####################################################################################################
####################       GENERATE UPSTREAM NETWORK      ##########################################
####################################################################################################
def getUpstreamDictMP(allUp, regionCOMIDs, PLUS):
    """
    Method: getUpstreamMP()
    Purpose: Chunk the comids in the region and find their upstream catchments on multiple cores.
    Params: allUp - dictionary of direct upstream links
            regionCOMIDs - list of COMIDs in the region
            PLUS - df of routing db
    Returns: fin - list of COMIDs that drain to catchments in region
    """
    #pull cores
    numCores = mp.cpu_count() - 2

    #chunk up data
    batch_size = int((len(regionCOMIDs) / numCores) + 1)
    chunk_iterator = []
    for i in range(numCores):
        mn, mx = i * batch_size, (i + 1) * batch_size
        args = regionCOMIDs[mn:mx], allUp, PLUS
        chunk_iterator.append(args)
    
    #call intersect function
    pool = mp.Pool(processes=numCores)
    results = pool.map(allUpstreamDict, chunk_iterator)
    pool.close()

    #concat results
    fin = {}
    for result in results: #list of dictionaries
        fin.update(result)
    return fin

def allUpstreamDict(args):
    """
    Method: allUpstream()
    Purpose: To find all upstream catchments in the list
    Params: args - tuple of:
                catchmnts - list of catchment IDs
                allUp - dictionary containing the direct upstream connections for each unique
                    FL_ComID in the routing database
                PLUS - dataframe of the routing database
    Returns: toSum - list of complete upstream connections (catchments and flowlines)
    """
    catchmnts, allUp, PLUS = args
    allCatch = {}
    for s in catchmnts:
        toSum = [int(s)] 
        new = [int(s)]
        while len(new) > 0:
            new, x = iterCheck(new, allUp, PLUS) 
            if len(x) > 0:
                toSum = list(set(toSum) - set(x))
                new = list(set(new) - set(x))
            new = list(set(new) - set(toSum)) #remove any items from new that are already stored in the network
            toSum = toSum + new
        allCatch[s] = list(set(toSum))
    return allCatch

def generateUpstream(PLUS):
    """
    Method: generateUpstream()
    Purpose: To compute the complete upstream network for the specified catchment
    Params:  PLUS - dataframe of the routing database
    Returns: toSum - list of complete upstream connections (catchments and flowlines)
    """
    print("Find direct upstream links...")
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

def iterCheck(coms, allUp, PLUS):
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
Global variables needing user input
VPUzone: NHD zone that the catchments are in (02 is MidAtlantic)
catchID: Name of column with unique catchment ID in the catchments shapefile
"""
VPUzone = "02" 
catchID = 'COMID'

"""
Method: main
Purpose: Call functions
"""
if __name__ == "__main__":
    #File Paths
    NHD_SHP = r"/media/imagery/HWGIT/NHDv2_CHWAv2.shp" #path to all catchments 
    regionCatch = r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/AccumulateNHDPlus/MD_Catchments/md_catchments.shp" #path to catchments in region
    sassPath = r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/ROUTING_DATA/SassDb_02_040721.csv" #path to routing DB
    NHD_lines = r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/ROUTING_DATA/CHWA_Lines.shp" #only needs if passing sas7bdat file for routing
    out_file = r"/media/imagery/usgs_sc/smcdonald/NHDv2_Agg/AccumulateNHDPlus/MD_Catchments/md_catch_paths_040821.shp" #path and name of output file

    writePath = True #true if you want to to add drainage path field; false if you just want to write all catchments draining to area

    #Read catchments into gdf
    print("Reading in catchments...")
    catchments = gpd.read_file(NHD_SHP)[[catchID, 'geometry']]
    #get list of COMIDs in region
    print("Reading catchments in region...")
    regionCOMIDS = gpd.read_file(regionCatch)
    regionCOMIDS = list(regionCOMIDS[catchID])
    #read in routing db
    allUp = {}
    PLUS = readRouting(sassPath, VPUzone, NHD_lines)
    #Generate direct upstream network
    allUp = generateUpstream(PLUS)
    #Determine all upstream catchments of region
    if writePath: # find and id unique flow paths
        print("Finding catchments draining to catchments in region...")
        resDict = getUpstreamDictMP(allUp, regionCOMIDS, PLUS)
        uniqueRegPaths = [] #find reg catchments that are not in drainage path of another reg catchment (duplicate paths)
        print("Finding unique drainage paths...")
        for reg in regionCOMIDS:
            exclude = False
            for d in resDict:
                if d != reg:
                    if reg in resDict[d]:
                        exclude = True
                        break
            if not exclude:
                uniqueRegPaths.append(reg)
        allCatch = []
        for u in uniqueRegPaths:
            allCatch += resDict[u]
        catchments = catchments[catchments[catchID].isin(allCatch)]
        catchments['drngPath'] = 0
        catchments['commonPath'] = None
        ct = 1
        for u in uniqueRegPaths:
            if len(resDict[u]) > 1: #if no drainage path (is sink) - leave 0
                if len(catchments[(catchments[catchID].isin(resDict[u])) & (catchments['drngPath'] != 0)]) > 0:
                    otherDrnPaths = list(set(list(catchments[(catchments[catchID].isin(resDict[u])) & (catchments['drngPath'] != 0)]['drngPath'])))
                    commonCatchments = list(catchments[(catchments[catchID].isin(resDict[u])) & (catchments['drngPath'] != 0)][catchID])
                    catchments.loc[catchments[catchID].isin(commonCatchments), 'drngPath'] = ct
                    ct += 1
                    cp = ''
                    for i in range(len(otherDrnPaths)):
                        cp += str(otherDrnPaths[i])
                        if i < len(otherDrnPaths) - 1:
                            cp += ','
                    catchments.loc[catchments[catchID].isin(commonCatchments), 'commonPath'] = cp
                    resDict[u] = list(set(resDict[u]) - set(commonCatchments)) #remove catchments that were just classed
                if len(resDict[u]) > 0:
                    catchments.loc[catchments[catchID].isin(resDict[u]), 'drngPath'] = ct
                    ct += 1
        print("Writing results...")
        catchments.to_file(out_file)
    else: # write out all catchments draining to area without path info
        print("Finding catchments draining to catchments in region...")
        allComsInReg = getUpstreamMP(allUp, regionCOMIDS, PLUS)
        #keep catchments in list
        catchments = catchments[catchments[catchID].isin(allComsInReg)]
        #write out results
        print("Writing results...")
        catchments.to_file(out_file)

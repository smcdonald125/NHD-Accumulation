import multiprocessing as mp
from shapely.geometry import Polygon, MultiPolygon, box
import numpy as np
import rasterio as rio
from rasterio.mask import mask
from rasterio.features import rasterize
import pandas as pd
import json
import geopandas as gpd

def aggregateRas_mp(shapes, stat, rasPath, rasType, convFact, cellSize, crs):
    """
    Method: aggregateRas_mp()
    Purpose: Divide catchments up to be multi-processed through statistic calculations.
    Params: shapes - dict of catchment IDs and geometries
            stat - statistic type to be calculated for the catchments (only sued for continuous data)
            rasPath - path to the raster data
            rasType - 0 or 1 specifying if raster is continuous or classed
            convFact - use defined conversion factor to apply to raster area
            cellSize - raster cell sized; needed to convert pixel count to area for classed data
            crs - vector crs from geopandas
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
    """
    cpus_minus_1 = mp.cpu_count() - 2
    if cpus_minus_1 == 0:
        cpus_minus_1 = 1

    pool = mp.Pool(processes=cpus_minus_1)

    chunk_iterator = []
    for i in range(len(shapes)): # one chunk per geom
        gdf_args = list(shapes.keys())[i], list(shapes.values())[i], rasPath, rasType, stat, convFact, cellSize, crs
        chunk_iterator.append(gdf_args)

    print("Number of Catchments: ", len(shapes))

    results = pool.map(aggregateRas, chunk_iterator) #list of dictionaries
    pool.close()

    result = {} #create one dict from list of dicts
    total = 0
    for d in results:
        total += len(d)
        result.update(d)
        if len(result) != total:
            print("Error: ", len(result), total)
    del results

    return result

def aggregateRas(args):
    """
    Method: aggregateRas()
    Purpose: calculate raster statistics for each catchment, find the upstream network for each catchment and
            accumulate the raster statistics.
    Params: key - catchment ID
            geom - catchment geometry
            rasPath - path to the raster data
            rasType - 0 or 1 specifying if raster is continuous or classed
            stat - statistic type to be calculated for the catchments (only needed for continuous data)
            convFact - user defined conversion factor to apply to raster area
            cellSize - raster cell sized; needed to convert pixel count to area for classed data
            crs - geopandas crs of vector zones
    Returns: sumDict - dictionary of catchment statistics and accumulated statistics
                        Continuous data: list where the catchment ID is the key and the data is a list of catchment statistics
                        Classed data: dataframe of class value and pixel count
    """
    key, geom, rasPath, rasType, stat, convFact, cellSize, crs = args
    sumDict = {}
    # snap and mask the raster by the zone
    # convert fiona geometry dictionary to shapely polygon
    t = geom['type']
    geom = geom['coordinates']
    while type(geom[0]) == list:
        geom = geom[0]
    try:
        geom = Polygon(geom)
    except Exception as e:
        print('----------------------------')
        print(e)
        print(type)
    ary, noData, area = mask_snap_raster(geom, rasPath, crs)
    #Calculate statistic for each catchment
    if rasType == 0: #continuous raster
        val, flag = agg_cont(ary, stat, noData)
        if flag:
            sumDict[key] = [float(val), float(area * convFact)]
    elif rasType == 1: #classed raster
            df = agg_classed(ary, noData, cellSize, convFact)
            if len(df) > 0: # replace with .empty?
                sumDict[key] = df.copy()
                del df
    return sumDict 
         
def agg_classed(ary, noData, cellSize, convFact):
    """
    Method: agg_classed()
    Purpose: Mask the classed raster to the current catchment and calculate the cell count
            for each class. Store data in a dataframe.
    Params: ary - numpy array of masked snapped raster
            noData - raster noData value
            cellSize - raster cell sized; needed to convert pixel count to area for classed data
            convFact - user defined conversion factor to apply to area
    Returns: df - pandas dataframe of raster values and counts
    """
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

def agg_cont(ary, statType, noData):
    """
    Method: agg_cont()
    Purpose: Mask the continuous raster to the current catchment and calculate the specified statistic
            and raster area within the catchment. The statistic options are: MAX, MIN, MEAN, MEDIAN and SUM.
    Params: ary - numpy array of masked snapped raster
            statType - string denoting the statistic type
            noData - raster noData value 
    Returns: catchSum - catchment statistic
    """
    try:
        #set no data value to nan for float rasters and 0 for integer rasters
        if statType != "SUM" and (type(ary[0][0][0]) == np.float32 or type(ary[0][0][0]) == np.float64):
            ary[ary == noData] = np.nan #set nodata to nan 
        else:
            ary[ary == noData] = 0
        if statType == "MAX":
            catchSum = float(np.nanmax(ary)) 
        elif statType == "MIN":
            catchSum = float(np.nanmin(ary)) 
        elif statType == "MEAN":
            catchSum = float(np.nanmean(ary)) 
        elif statType == "MEDIAN":
            catchSum = float(np.nanmedian(ary)) 
        elif statType == "SUM":
            catchSum = float(np.sum(ary))
        return catchSum, True #returns value
    except:
        return -1, False

def mask_snap_raster(geom, snap_raster_file, crs):
    """
    Method: mask_snap_raster()
    Purpose: Rasterize and snap a given geometry to a raster to be summarized by the given geometry.
    Params: geom - shapely geometry
            snap_raster_file - path to a raster file to mask the geometry by
            crs - geopandas crs from initial vector file
    Returns: masked_snap - numpy array that has been masked by the snapped geometry
    Authors: Labeeb Ahmed, Geographer, U.S. Geological Survey, lahmed@usgs.goc
             Sarah McDonald, Geographer, U.S. Geological Survey, smcdonald@usgs.gov
    """
    # 1 - Create Shapely BBox objects
    xmin, ymin, xmax, ymax = geom.bounds
    bbox = box(xmin, ymin, xmax, ymax)

    # 2 - Create a gdf out of bbox objects and use it to mask the snap raster
    mask_gdf = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=crs)
    mask_coords = [json.loads(mask_gdf.to_json())['features'][0]['geometry']]

    # 3 - Read in the snap raster and mask it
    overlap = True
    with rio.open(snap_raster_file, 'r') as snap:
        new_meta = snap.meta.copy()
        nodata = snap.nodatavals[0]
        
        # mask and create new transform
        try:
            masked_snap, new_transform = mask(snap, shapes=mask_coords, crop=True)
        except:
            overlap = False

    if overlap:
        # 4 - Convert polygons to array. This is collection of features to be rasterized
        shapes = [(geom, 1)]
        
        # 5 - Rasterize the feature collection
        height, width = masked_snap.shape[1:] if len(masked_snap.shape) > 2 else masked_snap.shape
        zone_array = rasterize(shapes, 
                                out_shape=(height, width), 
                                fill=0, 
                                transform=new_transform, 
                                all_touched=False)
        
        # 6 - Mask the snap raster by the rasterized geometry
        masked_snap = np.where(zone_array == 1, masked_snap, nodata)

        # 7 - Get area of zone
        area = ((zone_array == 1).sum()) # number of pixels that are in zone

        # return resulting array
        return masked_snap, nodata, area
    else:
        return np.array([nodata]), nodata, 0
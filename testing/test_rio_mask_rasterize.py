"""
Automate testing of Rasterio mask function vs Rasterio Rasterize function.

1. Loop through polygon features and count cells that are not no data
    - do not use p6 rasters (0 as nodata causes issues)
2. Get unique cell counts for the rasterized zones produced by rasterio already
3. Format data into a csv and join with local catchment shapefile to visualize
"""
import rasterio as rio
from rasterio import mask
import pandas as pd
import geopandas as gpd
import numpy as np

if __name__=="__main__":
    vector = r"C:/Users/smcdonald/Documents/Data/NHDv2_CHWAv2.shp"
    raster = r"C:/Users/smcdonald/Documents/Data/NHDv2_CHWAv2_rio.tif"
    out_csv = r"C:/Users/smcdonald/Documents/Data/rio_ras_v_mask.csv"

    # read in vector catchments
    print("Reading vector...")
    gdf = gpd.read_file(vector)
    gdf = gdf[['COMID_int', 'geometry']]

    # convert to list of tuples
    shapes = [ (row['geometry'], row['COMID_int']) for idx, row in gdf.iterrows()]
    del gdf
    print("Created shapes...")

    data = pd.DataFrame(columns=['COMID_int', 'Mask_Count'])
    # loop through shapes and mask data
    with rio.open(raster) as src:
        # get mask counts
        for shape in shapes:
            ary, t = rio.mask.mask(src, [shape[0]], all_touched=False, invert=False, nodata=0, filled=True, crop=True)
            unique, counts = np.unique(ary, return_counts=True)
            unique = list(unique)
            counts = list(counts)
            if 0 in unique:
                idx = unique.index(0)
                del counts[idx]
            c = sum(counts)
            data.loc[len(data)] = [shape[1], c]
        
        print(data)
        # get raster zone counts
        ary = src.read(1)
        unique, counts = np.unique(ary, return_counts=True)
        del ary
        unique = list(unique)
        counts = list(counts)
        df = pd.DataFrame(data={'COMID_int':unique, 'Ras_Count':counts})
        print(df)
        data = data.merge(df, on='COMID_int', how='outer')
        del df
    
    # write results
    data.to_csv(out_csv, index=False)
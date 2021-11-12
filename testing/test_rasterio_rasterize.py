import rasterio as rio
from rasterio.features import rasterize
from rasterio.enums import Resampling
import geopandas as gpd
import math
import affine
import numpy as np

if __name__=="__main__":
    vector = r"C:/Users/smcdonald/Documents/Data/NHDv2_CHWAv2.shp"
    out_raster = r"C:/Users/smcdonald/Documents/Data/NHDv2_CHWAv2_rio.tif"
    template_raster = r"G:/ImageryServer/1m_LU_2017/version1/phase6_10m/IR_2017_10m.tif" # used for metadata
    cellSize = 10.0

    # read in vector catchments
    print("Reading vector...")
    gdf = gpd.read_file(vector)
    gdf = gdf[['COMID_int', 'geometry']]

    # get extent of raster
    bounds = list(gdf.total_bounds) # minx , miny , maxx , maxy
    sh = (math.ceil((bounds[3] - bounds[1])/cellSize), math.ceil((bounds[2] - bounds[0])/cellSize)) # height, width
    print("Shape: ", sh)

    # create transform from gdf
    transform = affine.Affine(cellSize, 0, bounds[0], 0, -cellSize, bounds[3])
    print("Transform: ", transform)

    # convert to list of tuples
    shapes = [ (row['geometry'], row['COMID_int']) for idx, row in gdf.iterrows()]
    del gdf
    print("Created shapes...")

    # convert polygons to an array
    ary = rasterize(shapes, out_shape=sh, fill=0,  transform=transform, all_touched=False)
    del shapes
    print("Created Array... ", ary.shape, np.amin(ary), np.amax(ary))

    # create metadata for raster
    with rio.open(template_raster) as src:
        meta = src.meta
    meta.update({'nodata':0,
                'dtype':'uint32',
                'width':sh[1],
                'height':sh[0],
                'transform':transform})

    # write out array to geotiff
    print("Writing to tiff...")
    with rio.open(out_raster, 'w', **meta, compress="LZW") as dataset:
        try:
            dataset.write(ary)
        except:
            try:
                dataset.write(ary, 1)
            except Exception as e:
                print("Write To Tiff Failed")
                print(e)
    del ary
    
    # create pyramids
    factors = [2, 4, 8, 16, 32, 64, 128, 256, 512]
    with rasterio.open(out_raster, 'r+') as dst:
        dst.build_overviews(factors, Resampling.nearest)
        dst.close()
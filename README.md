# NHD-Accumulation
Open-source Python code to accumulate raster data by NHD+ v2 catchments and by high res catchments.

The two scripts are standalone, open-source, python 3.7 scripts. Both share the same required libraries, listed below. The libraries
can be installed using a conda environment with the command: conda install *libraryName*
  
  o	pandas
  
  o gdal
  
  o	geopandas
  
  o	numpy
  
  o	rasterio
  
  o	shapely
  
  o	fiona 

## NHDPlusV2

This script requires: 

 * NHDPlusV2 catchments
 
 * NHDPlusV2 flowlines (needed to remove upstream coastline connections)
 
 * routing database, unzipped nhdplusv2_us_sas7bdat.zip file found here: https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47
    * The routing database can be read in its SAS format, or can be read in as a CSV. The database covers CONUS, and should be reduced to a specified zone if possible.
  The code will produce a CSV of the database with the specified VPU (Vector-Processing Unit) only, if the VPUzone field is set. The produced CSV will only contain
  information required to build the upstream networks.

## NHD-HR

This script can use either NHD High Res catchments and PlusFlow table, or ecoSHEDS catchments (optional truncated flowlines). 
The isECO field needs to be set to True if using ecoSHEDs and False if using NHD High Res.

## Accumulated Data

Both scripts can accumulate either a continuous raster or a classified raster. AccumulateNHDPlus can accumulate tabular data, as long as there is a field with the NHD COMID. 

### Continuous Raster

These data can be summarized using one of the following statistics:

  MIN,
  MAX,
  MEAN,
  MEDIAN,
  SUM
  
The result of the statistic, as well as "raster area" (total area according to the raster), will be accumulated for each catchment. 

### Classified Raster

These data will have the total pixel count for each unique class, excluding NoData, recorded and accumulated. The total accumulated
area for each class is calculated by mutiplying the pixel count with a user-defined conversion factor.

### Tabular Data

These data will have the accumulated values for user specified fields, with the final output having "Ws" in front of the original column name.

### Contact

Sarah McDonald, Geographer, US Geological Survey
Chesapeake Bay Program Office
Lower Mississippi-Gulf Water Science Center
smcdonald@chesapeakebay.net 

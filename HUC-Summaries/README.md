# Summarize Tabular Data by HUC
This script utilizes the routing information in the HUC data to build a HUC-scale network. Metrics to accumulate are
assumed to be in the same file as the vector HUC layer, with the routing field. Data can be accumulated up or 
down stream.

## How to Use
### Library Set Up
Base and open-source Python libraries are required, including pandas and geopandas. 
- os
- sys
- logging
- time
- pathlib
- pandas
- geopandas

### User Input
This script has 6 required user-input arguments. Unlimited optional arguments are also available.
1. Path to routing file.
2. Direction to build network, options are upstream or downstream
3. Routing type. Enter feature if the routing network uses the unique feature ids. Enter node if the network uses nodes to connect feature ids
4. Unique feature ID
5. The 'to field' in the routing data
6. The 'from field' in the routing data
7+. Fields to accumulate, other than area. List as many as desired, separated by spaces.

### Example Run
To summarize average bank height (avg_bank_height) by HUC12 (stored in HUC12.gpkg) downstream, where the routing information is by unique huc12
in the field tohuc.
- python accumulateHUCs.py HUC12.gpkg downstream huc12 tohuc huc12 avg_bank_height

## Output Files
This script generates 3 output files, all stored in the same directory as the input file:
1. GPKG of summarized data appended to original data. File is same name as original, with "_upstream" or "_downstream" appended.
2. Comma-separated file of full upstream/downstream network per HUC. File is same name as output GPKG with "_network".
3. Log file of the script run. File is named accumulateHUCs_log.txt.

# Contact
Author: Sarah M. McDonald, Geographer
Affiliations: U.S. Geological Survey, Lower Mississippi-Gulf Water Science Center. Chesapeake Bay Program.
Contact: smcdonald@chesapeakebay.net

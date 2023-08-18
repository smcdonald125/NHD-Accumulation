"""
Script: accumulateHUCs.py
Purpose: This script uses routing information to build upstream and downstream networks and accumulates
            area in acres, the number of features upstream/downstream, and user-specified metrics. 
Requirements: The code expects vector layer with geometries and routing information.
User arguments: The script expects 6 command line arguments from the user, with an optional 7th:
        1. Path to routing file.
        2. Direction to build network, options are upstream or downstream
        3. Routing type. Enter feature if the routing network uses the unique feature ids. Enter node if the network uses nodes to connect feature ids
        4. Unique feature ID
        5. The 'to field' in the routing data
        6. The 'from field' in the routing data
        7. Fields to accumulate, other than area. List as many as desired, separated by spaces.
Example: To run a downstream accumulation for HUC12s where the script and data are in the same directory:
            python accumulateHUCs.py HUC12.gpkg downstream feature huc12 tohuc huc12

        - The results will exist in the same folder as the HUC12.gpkg with the original name plus direction. For this example,
          the resulting gpkg will be HUC12_downstream.gpkg.
        - the network will be stored in the same directory, with the original name plus direction. For this example,
          the resulting network file is comma-separated and named HUC12_downstream_network
        - the script logs the code run in the same directory under the name accumulateHUCs_log
Author: Sarah McDonald, Geographer, U.S. Geological Survey, Chesapeake Bay Program
Contact: smcdonald@chesapeakebay.net
"""

import os
import sys
from pathlib import Path
import logging
import time
import pandas as pd 
import geopandas as gpd 

logger = logging.getLogger(__name__)

def read_commandLine() -> tuple:
    """
    read_commandLine This method reads in user input and validates them.

    Returns
    -------
    tuple
        routing_path : str
            string path to routing file.
        direction : str
            direction to build network (upstream or downstream).
        routing_type : str
            feature if the network uses feature IDs to build the network
            node if the network uses nodes to build connection by feature IDs
        uniqueID : str
            field name of the unique feature IDs
        toFld : str
            field name that is the "to" in the routing
        fromFld : str
            field name that is the "from" in the routing
        metrics : str
            list of metrics to accumulate
    """
    # read command line arguments
    args = sys.argv
    if len(args) >= 7:
        # user input
        routing_path = args[1]
        direction = args[2]
        routing_type = args[3]
        uniqueID = args[4]
        toFld = args[5]
        fromFld = args[6]
        metrics = []
        if len(args) > 7:
            for i in range(7, len(args)):
                metrics.append(args[i])

        # set up logger
        logging_setup(Path(routing_path).parent)

        # read routing
        if not os.path.isfile(routing_path):
            logger.warning(f"Routing Path does not exist: {routing_path}")
            sys.exit()
        else:
            header = list(gpd.read_file(routing_path, rows=0))

        if direction not in ['upstream', 'downstream']:
            logger.warning(f"Expected second argument (direction) to be upstream or downstream: {direction}")
            sys.exit()
        if routing_type not in ['feature', 'node']:
            logger.warning(f"Expected third argument (routing_type) to be feature or node: {routing_type}")
            sys.exit()
        if uniqueID not in header:
            logger.warning(f"Fourth argument (uniqueID) does not exist in the routing file {uniqueID}\nOptions are {header}")
            sys.exit()
        if toFld not in header:
            logger.warning(f"Fifth argument (to field) does not exist in the routing file {toFld}\nOptions are {header}")
            sys.exit()
        if fromFld not in header:
            logger.warning(f"Sixth argument (from field) does not exist in the routing file {fromFld}\nOptions are {header}")
            sys.exit()
        if len( list(set(metrics) - set(header)) ) > 0:
            missing = list(set(metrics) - set(header))
            logger.warning(f"Seventh argument (metrics to accumulate) does not exist in the routing file {missing}\nOptions are {header}")
            sys.exit()
    else:
        m = f"Missing arguments. Expected:\n\t1. Path to routing file."
        m += f"\n\t2. Direction to build network, options are upstream or downstream"
        m += f"\n\t3. Routing type. Enter feature if the routing network uses the unique feature ids. Enter node if the network uses nodes to connect feature ids"
        m += f"\n\t4. Unique feature ID"
        m += f"\n\t5. The 'to field' in the routing data"
        m += f"\n\t6. The 'from field' in the routing data"
        m += f"\n\t7. Optional: metrics to accumulate (area is by default)"
        logging.error(m)
        sys.exit()

    logger.info(f"routing_path: {routing_path}")
    logger.info(f"direction: {direction}")
    logger.info(f"routing_type: {routing_type}")
    logger.info(f"uniqueID: {uniqueID}")
    logger.info(f"toFld: {toFld}")
    logger.info(f"fromFld: {fromFld}")
    logger.info(f"metrics: {metrics}")

    return direction, routing_type, uniqueID, toFld, fromFld, routing_path, metrics

def logging_setup(folder):
    # initiate logging
    logging.captureWarnings(True)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s %(message)s',
                        handlers=[logging.FileHandler(f"{folder}/accumulateHUCs_log.txt", mode='a'),
                                  stream_handler])
    
    logger.info(f"\n------------------------------------------\n------------------------------------------\n")

def generateDirectLinks_node(routingDB, direction, toFld, fromFld, uniqueID):
    """
    Method: generateDirectLinks_node()
    Purpose: To compute the complete network for the specified feature.
    Params:  routingDB - dataframe of the routing database
             direction - upstream or downstream
             toFld - to node field
             fromFld - from node field
             uniqueID - unique feature ID
    Returns: allUp - dictionary of direct links
    """
    # TODO - verify node data don't have same issue as feature (to or from field is same as unique ID)
    # Create dict of direct upstream connections -- faster than referencing table directly each time
    allUp = {}
    for idx, row in routingDB.iterrows():
        if direction == 'downstream':
            fromN = row[fromFld]
            x = list(routingDB[routingDB[toFld] == fromN][uniqueID])[0]
        else:
            fromN = row[toFld]
            x = list(routingDB[routingDB[fromFld] == fromN][uniqueID])
        allUp[row[uniqueID]] = x.copy()
    return allUp

def generateDirectLinks_feature(routingDB, direction, toFld, fromFld, uniqueID):
    """
    Method: generateDirectLinks_feature()
    Purpose: To compute the complete network for the specified feature.
    Params:  routingDB - dataframe of the routing database
             direction - upstream or downstream
             toFld - to node field
             fromFld - from node field
             uniqueID - unique feature ID
    Returns: allUp - dictionary of direct links
    """
    # Create dict of direct upstream connections -- faster than referencing table directly each time
    allUp = {}
    if direction == 'downstream':
        cnnField = fromFld
    else:
        cnnField = toFld

    if cnnField == uniqueID: # if IDs are the same, loop through records and use opposite field
        if direction == 'downstream':
            cnnField = toFld
        else:
            cnnField = fromFld
        for idx, row in routingDB.iterrows():
            if row[cnnField] in allUp: # add record to list
                allUp[row[cnnField]] += [row[uniqueID]]
            else: # create record
                allUp[row[cnnField]] = [row[uniqueID]]
    else:
        for idx, row in routingDB.iterrows():
            allUp[row[uniqueID]] = [row[cnnField]]

    return allUp

def generateNetwork(allUp, s):
    """
    Method: generateNetwork()
    Purpose: To compute the complete upstream network for the specified catchment
    Params: allUp - dictionary containing the direct upstream connections for each unique
                    FL_ComID in the routing database
            s - the current catchment ID
            routingDB - dataframe of the routing database
    Returns: toSum - list of complete upstream connections (catchments and flowlines)
    """
    network = [s] 
    if s in allUp: # skip if no direct links
        new = [s]
        while len(new) > 0:
            new, x, allUp = iterCheck(new, allUp) 
            if len(x) > 0:
                network = list(set(network) - set(x))
                new = list(set(new) - set(x))
            new = list(set(new) - set(network)) #remove any items from new that are already stored in the network
            network = network + new
    return list(set(network))

def iterCheck(coms, dirLinks):
    """
    Method: iterCheck()
    Purpose: find direct upstream connections for all FL_ComIDs in a list and build a list
            of the unique IDs found. Called by generateNetwork
    Params: coms - list of IDs to find direct connections of
            dirLinks - dictionary containing the direct connections for each unique
                    ID in the routing database
    Returns: new - list of unique upstream catchment IDs
            notFound - list of unique IDs that did not exist in allUp
    """
    new, notFound = [], []
    for c in coms:
        if c in dirLinks:
            if len(dirLinks[c]) > 0:
                new = new + dirLinks[c]
        else:
            """ 
            This warning can have several causes. Verify the feature is truly not connected for one of the reasons below.
                1. Insufficient routing information (is the full network present in the routing data?)
                2. The feature is a headwater (downstream accumulation) or an outlet (upstream accumulation)
                3. The feature is a sink - no connections in or out?
            """
            logger.warning(f"\t\t{c} - no connections")
            notFound.append(c)
            dirLinks[c] = []
    return list(set(new)), notFound, dirLinks

def accumulate(directLinks, routingDB, uniqueID, network, metrics):
    """
    accumulate builds full upstream/downstream network and sums acres and calculated number of features in network.

    Parameters
    ----------
    directLinks : dict
        dictionary of direct connections
    routingDB : gpd.GeoDataFrame
        database of routing and data
    uniqueID : str
        unique feature ID
    network : dict
        dictionary of network; empty dict if network needs to be generated
    metrics : list
        metrics to accumulate

    Returns
    -------
    tuple
        network : dict
            full upstream or downstream network
        routingDB : gpd.GeoDataFrame
            data with 2 new columns, acresWS and networkCnt. acresWs is the sum of acres field and networkCnt is 
            number of features in the network.
    """
    # dictionary of full network
    network = {}

    # set unique ID to index
    routingDB = routingDB.set_index(uniqueID)

    # create list of accumulated metric field names
    metricsWs = [f"{m}Ws" for m in metrics]
    
    # loop through IDs
    for i in routingDB.index:

        if i not in network:
            # build network
            network[i] = generateNetwork(directLinks, i)

        # sum values
        routingDB.loc[i, metricsWs] = list(routingDB[routingDB.index.isin(network[i])][metrics].sum())
        routingDB.loc[i, 'networkCnt'] = len(network[i])

    # return data
    return network, routingDB


if __name__=="__main__":
    # retrieve user input
    logger.info("Validating input")
    user_input = read_commandLine()
    direction, routing_type, uniqueID, toFld, fromFld, routing_path, metrics = user_input

    # building output paths
    p = Path(routing_path)
    outputDB = f"{p.parent}/{p.name.split('.')[0]}_{direction}.gpkg"
    network_path = f"{p.parent}/{p.name.split('.')[0]}_{direction}_network"

    # check if output exists and prompt to overwrite
    if os.path.isfile(outputDB):
        val = input(f"\n\nOutput data exists {outputDB}\nDo you want to overwrite (y or n)?: ")
        while val not in ['y', 'n']:
            val = input(f"Please enter y for yes or n for no: ")

        logger.info(f"Output data exists {outputDB}")
        if val == 'n':
            logger.info(f"User selected to not overwrite. Exiting.")
            sys.exit(0)
        else:
            logger.info(f"User selected to overwrite. Continuing.")

    # prep data
    logger.info("Reading data and calculating area in acres")
    routingDB = gpd.read_file(routing_path)
    logger.info(f"There are {len(routingDB)} records")
    routingDB.loc[:, 'acres'] = routingDB.geometry.area / 4046.86 # TODO: assuming square meters
    metrics.append('acres')

    # check for existing network
    network, directLinks = {}, {}
    if os.path.isfile(network_path):

        logger.info(f"Network file exists. Reading in {network_path}")
        network_df = pd.read_csv(network_path)
        network_df = network_df.set_index(list(network_df)[0]) # first column is key
        for idx, row in network_df.iterrows():
            network[idx] = list(row)
        del network_df

    else:

        # build direct upstream or downstream
        logger.info("Generating direction connections")
        if routing_type == 'node':
            directLinks = generateDirectLinks_node(routingDB, direction, toFld, fromFld, uniqueID)
        else:
            directLinks = generateDirectLinks_feature(routingDB, direction, toFld, fromFld, uniqueID)

    # accumulate data
    logger.info(f"Accumulating {metrics}")
    network, routingDB = accumulate(directLinks, routingDB, uniqueID, network, metrics)
    routingDB.loc[:, 'networkCnt'] = routingDB.networkCnt.astype(int)

    # write data
    logger.info(f"Writing results to {outputDB}")
    routingDB.to_file(outputDB, layer='main', driver="GPKG")

    # write network
    if not os.path.isfile(network_path):
        logger.info(f"Writing network to {network_path}")
        
        # need max length of network to create header, otherwise may cause read errors
        max_value = 0
        for i in network:
            if len(network[i]) > max_value:
                max_value = len(network[i])

        with open(network_path, "w") as f:
            # create header to match longest list
            header = "0"
            for i in range(1, max_value+1):
                header = f"{header},{i}"
            f.write(header)

            # add line for each feature's network
            for n in network:
                m = f"{n}"
                for u in network[n]:
                    m = f"{m},{u}"
                m = f"{m}\n"
                f.write(m)

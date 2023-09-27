##################################################################################
# Study: Street network generation
# Purpose: subdivide roads information by cells for parallelization of road width calculus
# Author: Marjolaine Lannes
# Creation date: January 17, 2023
# Note: Get roads information (id, geometry) of OSM roads for each cell of the grid
##################################################################################
import pandas as pd
from shapely.geometry import Polygon, Point, LineString
import pickle, sys
import geopandas as gpd
import time
from itertools import repeat
import numpy as np
from shapely.wkt import loads

# Input/output files
path = "//"
roads_f = path + "temp/load_OSM/roads_links.csv"
road_ids_file = path + "temp/load_OSM/road_cells_IDF.dat"
roads_output_file = path + "temp/load_OSM/roads_in_cells.csv"

# Load data
roads_df = pd.read_csv(roads_f)
roads_df.astype({'GEOMETRY':'str'})
roads_df.GEOMETRY = roads_df.GEOMETRY.apply(loads)
roads = gpd.GeoDataFrame(roads_df,crs='epsg:2154', geometry = 'GEOMETRY')
n_roads = roads.shape[0]
print("Number of roads: ",n_roads)

# Initialize road_cells Dataframe for the next algorithm
road_ids = {'ID':[], 'CELLS':[]}
road_ids['ID'] = roads_df['ID']
road_ids['CELLS'] = [[] for i in repeat(None, len(road_ids))]
print(roads_f + ' loaded')

# Save data for the roads within each cell
start = time.time()
# Load grid file
grid_file = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
outdir = path + 'temp/load_OSM/subdomains_roads_lamb93/'
grid = pickle.load(open(grid_file, 'rb'))

for poly_num, poly in enumerate(grid['polygon']):
    timer_poly = time.time()
    # Get the cell information
    x = grid['we_id'][poly_num]
    y = grid['sn_id'][poly_num]
    ## one file per cell with all roads included in it
    prefix = 'we_id_from_' + str(x) + '_sn_id_from_' + str(y)
    foutname = outdir + 'sub_' + prefix + '.dat'
    # Select roads intersecting the cell
    subdata = roads[roads.geometry.intersects(poly)]
    print('Road data croped in cell...')
    # For each road in the cell, get roads data and save cell number in the network file
    roads_in_cell = {'ROAD_ID': list(subdata.ID), 'GEOMETRY': list(subdata.GEOMETRY), 'LINE': list(subdata.LINE),
                     'NODE_A': list(subdata['NODE_A']), 'NODE_B': list(subdata['NODE_B']), 'LENGTH':list(subdata.LENGTH),
                     'XA': list(subdata.XA), 'XB': list(subdata.XB), 'YA': list(subdata.YA), 'YB': list(subdata.YB)}
    for road_id in roads_in_cell['ROAD_ID'] :
        r_index = roads_df[roads_df.ID == road_id].index[0] # road index in network file
        cells_list = road_ids['CELLS'][r_index] # get existing list of cells for this road
        cells_list.append(poly_num)
        road_ids['CELLS'][r_index] = cells_list
    # Save it
    print('Now saving....')
    pickle.dump(roads_in_cell, open(foutname, 'wb'))
    #Print information
    end = time.time()
    print('Done with cell ', poly_num, "in ", (end - timer_poly), " seconds" )
    if poly_num % 100 == 0 :
        print("time for cell ", poly_num, " roads:", (end - start), " seconds" )

# Save the lists of cells for each road
pickle.dump(road_ids, open(road_ids_file, 'wb'))
network_cells = pd.DataFrame(road_ids)
df = pd.merge(left=roads, right=network_cells, how='left', on=['ID'])
df.to_csv(roads_output_file)
print('Cells added to road data...')


##################################################################################
# Study: Street network generation
# Purpose: to save roads geometry from OSM
# Author: Marjolaine Lannes
# Creation date: January 19, 2023
# Note: Save the distribution of each road between the cells
# Output (*) columns: ['ROAD_ID', 'NODE_A', 'NODE_B', 'XA', 'YA', 'XB', 'YB', 'LENGTH', 'ROUNDABOUT', 'TUNNEL', 'LINE', 'GEOMETRY', 'CELLS', 'DISTRIBUTION']
##################################################################################
import pickle
import pandas as pd
import geopandas as gpd
import time
from shapely.wkt import loads

# Input/output files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
roads_file = path + "temp/data/load_OSM/roads_in_cells.dat"
road_cells_file = path + "temp/data/load_OSM/road_cells_IDF.dat"
roads_per_cell_f = path + "temp/data/load_OSM/roads_list_per_cell.dat"

# Load roads data
grid=pickle.load(open(gridfile,'rb'))
#road_cells = pickle.load(open(road_cells_file,'rb'))

# For each OSM road, calculate the percentage attributed to each cell it intersects
roads_in_cell = {}
#for poly_num in range(8250):
for poly_num in [4672]:
    start = time.time()
    # Select data in the cell
    x = grid['we_id'][poly_num]
    y = grid['sn_id'][poly_num]
    prefix = 'we_id_from_' + str(x) + '_sn_id_from_' + str(y)
    cell_roads_f = path + 'temp/load_OSM/subdomains_roads_lamb93/sub_' + prefix + '.dat'
    cell_roads_data = pickle.load(open(cell_roads_f,'rb'))
    cell_roads_df = pd.DataFrame(cell_roads_data)
    cell_roads_df = cell_roads_df.astype({'ROAD_ID':str})
    roads_in_cell[poly_num] = list(cell_roads_df['ROAD_ID'])
    if poly_num == 4672 :
        print(len(roads_in_cell[poly_num]))
    end = time.time()
    print("Cell number ", poly_num, " in ", end-start," seconds")

# Save it as .dat and GeoDataFrame
print("Calculations done ! Now saving...")
pickle.dump(roads_in_cell, open(roads_per_cell_f, 'wb'))

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
roads_output_dict_file = path + "temp/data/load_OSM/roads_cells_distrib_dict.dat"
roads_per_cell_f = path + "temp/data/load_OSM/roads_list_per_cell.dat"

# Load roads data
roads_df = pd.read_csv(roads_file)
grid=pickle.load(open(gridfile,'rb'))
road_cells = pickle.load(open(road_cells_file,'rb'))
road_cells_distribution = {}


# For each OSM road, calculate the percentage attributed to each cell it intersects
roads_df.astype({'GEOMETRY':'str'})
roads_df.GEOMETRY = roads_df.GEOMETRY.apply(loads)
roads = gpd.GeoDataFrame(roads_df,crs='epsg:2154', geometry = 'GEOMETRY')
roads_in_cell = {}
for poly_num, poly in enumerate(grid['polygon']):
    start = time.time()
    # Select data in the cell
    x = grid['we_id'][poly_num]
    y = grid['sn_id'][poly_num]
    prefix = 'we_id_from_' + str(x) + '_sn_id_from_' + str(y)
    cell_roads_f = path + 'temp/data/load_OSM/subdomains_roads_lamb93/sub_' + prefix + '.dat'
    cell_roads_data = pickle.load(open(cell_roads_f,'rb'))
    cell_roads_df = pd.DataFrame(cell_roads_data)
    roads_in_cell[poly_num] = list(cell_roads_df['ROAD_ID'])
    print('Data loaded...')
    # Calculate the proportion of the street in the cell for each street.
    for road_id in cell_roads_df['ROAD_ID'] :
        road_index = roads_df[roads_df.ID == road_id].index[0]
        cells_list = road_cells['CELLS'][road_index]
        if len(cells_list) == 1 :
            road_cells_distribution['DISTRIBUTION'][road_index][0] = 1
        else :
            geom = roads.loc[road_index,'GEOMETRY']
            line_length = roads_df.loc[road_index,'LENGTH']
            inter = poly.intersection(geom).length
            distrib = inter / line_length
            i = cells_list.index(poly_num)
            road_cells_distribution['DISTRIBUTION'][road_index][i] = distrib
    end = time.time()
    print("Cell number ", poly_num, " in ", end-start," seconds")

# Save it as .dat and GeoDataFrame
print("Calculations done ! Now saving...")
pickle.dump(road_cells_distribution, open(roads_output_dict_file, 'wb'))
pickle.dump(roads_in_cell, open(roads_per_cell_f, 'wb'))

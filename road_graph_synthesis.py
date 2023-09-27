##################################################################################
# Study: Street network generation
# Purpose: Synthesize roads information (including their width)
# Author: Marjolaine Lannes
# Creation date: February 28
# Note: Get width data from each cell and calculate each road weighted mean width
# Output (*) columns: ['ROAD_ID', 'NODE_A', 'NODE_B', 'XA', 'YA', 'XB', 'YB', 'LENGTH', 'ROUNDABOUT', 'TUNNEL', 'LINE', 'GEOMETRY', 'CELLS', 'DISTRIBUTION', 'WIDTH']
##################################################################################
import pickle
import pandas as pd
import numpy as np

# Input files
path = "//"
f_roads_distribution = path + "temp/load_OSM/roads_in_cells_distribution.dat"
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
road_widths_dir = path + 'temp/roads_graph/subdomains_roads_widths/'
foutput = path + 'temp/roads_graph/road_attributes.dat'

# Load data
#roads_data = pickle.load(open(f_roads_distribution,'rb'))
roads_df = pd.read_csv(f_roads_distribution)
#roads_df = pd.DataFrame(roads_data)
grid=pickle.load(open(gridfile,'rb'))
print('Data loaded...')

# nodes_OSM['LINKS'] = nodes_OSM['LINKS'].apply(eval) : tester ici ???

# Initialize output data
road_properties_detailled = {'ROAD_ID': [], 'CELLS':[],'WIDTHS':[], 'WIDTH':[]}
road_properties_detailled['ROAD_ID'] = roads_data['ROAD_ID']
n_roads = len(road_properties_detailled['ROAD_ID'])
for i in range(n_roads) :
    read_cells_list = roads_data['CELLS'][i]
    if type(read_cells_list) == str:
        read_cells_list = read_cells_list[1:-1].split(', ')
        road_cells = [int(x) for x in read_cells_list]
    else :
        road_cells = [int(x) for x in read_cells_list]
    road_properties_detailled['CELLS'].append(road_cells)
    n_cells = len(road_cells)
    properties_list = [0] * n_cells
    road_properties_detailled['WIDTHS'].append(properties_list)
print("Width initialized")
# road_properties_detailled['CELLS'] = roads_data['CELLS']
# for road_cells in roads_data['CELLS'] :
#     n_cells = len(road_cells)
#     properties_list = [0] * n_cells
#     road_properties_detailled['WIDTHS'].append(properties_list)

# For each poly, get their roads' attributes
for poly_num, poly in enumerate(grid['polygon']) :
    # input file for the cell
    x0 = grid['we_id'][poly_num]
    y0 = grid['sn_id'][poly_num]
    prefix = 'sub_we_id_from_' + str(x0) + '_sn_id_from_' + str(y0) + '.dat'
    # load roads widths for the cell
    road_widths_data = pickle.load(open(road_widths_dir + prefix, 'rb'))
    road_widths_df = pd.DataFrame(road_widths_data)
    # save road properties (one list of widths per road with all data from cells intersected by this road)
    for road_index_width, road_id in enumerate(road_widths_df['ROAD_ID']):
        road_width = road_widths_df.loc[road_index_width, "mean_width"]
        road_index = roads_df[roads_df.ROAD_ID == road_id].index[0]
        cells_list = road_properties_detailled['CELLS'][road_index]
        # if type(cells_list) == str :
        #     cells_list = cells_list[1:-1].split(', ')
        #     cell_index = cells_list.index(str(poly_num))
        # else :
        #     cell_index = cells_list.index(poly_num)
        cell_index = cells_list.index(poly_num)
        road_properties_detailled['WIDTHS'][road_index][cell_index] = road_width
print("Width calculated")

# Then, average properties for roads intersecting multiple cells
for road_index, road in enumerate(roads_data['ROAD_ID']):
    # polys = roads_data['CELLS'][road_index]
    # distribution = roads_data['DISTRIBUTION'][road_index]
    read_distribution = roads_data['DISTRIBUTION'][road_index]
    if type(read_distribution) == str:
        read_distribution_list = read_distribution[1:-1].split(', ')
        distribution = [int(x) for x in read_distribution]
    else :
        distribution = [int(x) for x in read_distribution]
    widths = road_properties_detailled['WIDTHS'][road_index]
    if len(widths) != len(distribution) :
        print(road, widths, distribution, read_distribution)
    road_properties_detailled['WIDTH'] = np.average(widths, weights=distribution)
print("Width averaged")

# Save road attributes
road_attributes = pd.DataFrame.from_dict(road_properties_detailled)
road_attributes = road_attributes.drop(columns=['CELLS','WIDTHS'])
road_attributes = pd.merge(left=roads_df, right=road_attributes, how='left', on=['ROAD_ID']) # or 'ROAD_ID ?
pickle.dump(road_attributes, open(foutput,'wb'))
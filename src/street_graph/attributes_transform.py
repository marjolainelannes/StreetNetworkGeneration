##################################################################################
# Study: Street network generation
# Purpose: Attribute properties to streets in the Île-de-France graph
# Note: Get height/width data from each cell and calculate each street weighted mean height/width
##################################################################################
import pickle, os, time
import pandas as pd
import numpy as np
from ast import literal_eval

# Parameters
save_streets_on_grid = False
save_canyon_attribute = True

# Input files
path = "../../"
temp_dir = path + "temp/"
data_dir = path + "data/"
out_dir = path + "output/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Files
f_streets = temp_dir + "street_graph/graph_transform/street_links_distrib.csv"
gridfile = data_dir + 'grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
street_height_dir = temp_dir + 'street_graph/attributes/subdomains_streets_height_width/'
foutput = out_dir + 'street_attributes'

# Load data
init = time.time()
streets_df = pd.read_csv(f_streets)
streets_df['CELLS'] = streets_df['CELLS'].apply(lambda x: literal_eval(x))
streets_df['DISTRIBUTION'] = streets_df['DISTRIBUTION'].apply(lambda x: literal_eval(x))
n_streets = streets_df.shape[0]
grid = pickle.load(open(gridfile,'rb'))
print('Data loaded...')

# Initialize data
street_properties_detailled = {'STREET_ID': [], 'CELLS':[], 'DISTRIBUTION':[], 'N_CELLS':[] ,'WIDTHS':[], 'LEFT_HEIGHTS':[], 'RIGHT_HEIGHTS':[]}
street_properties_detailled['STREET_ID'] = streets_df['STREET_ID']
street_properties_detailled['CELLS'] = streets_df['CELLS']
street_properties_detailled['DISTRIBUTION'] = streets_df['DISTRIBUTION']
for street_cells in streets_df['CELLS'] :
    n_cells = len(street_cells)
    properties_list = [0] * n_cells
    street_properties_detailled['N_CELLS'].append(n_cells)
    street_properties_detailled['WIDTHS'].append(properties_list.copy())
    street_properties_detailled['LEFT_HEIGHTS'].append(properties_list.copy())
    street_properties_detailled['RIGHT_HEIGHTS'].append(properties_list.copy())
time_end = time.gmtime(time.time() - init)
line = str(time_end.tm_mday) + "days"+ str(time_end.tm_hour) + 'hour' + str(time_end.tm_min) + 'min' + str(time_end.tm_sec)+ 'sec.'
print("Data initialized - time", line)
    
# For each poly, get their streets' attributes
for poly_num, poly in enumerate(grid['polygon']) :
    # input files
    x0 = grid['we_id'][poly_num]
    y0 = grid['sn_id'][poly_num]
    prefix = 'sub_we_id_from_' + str(x0) + '_sn_id_from_' + str(y0) + '.dat'
    # load data in the cell
    street_heights_data = pickle.load(open(street_height_dir + prefix, 'rb'))
    street_heights_df = pd.DataFrame(street_heights_data)
    # Get street properties
    for street_SUB_index, street_id in enumerate(street_heights_data['STREET_ID']):
        # Street width/height
        street_width = street_heights_df.loc[street_SUB_index, "WIDTH"]
        street_height_left = street_heights_df.loc[street_SUB_index, "HEIGHT_LEFT"]
        street_height_right = street_heights_df.loc[street_SUB_index, "HEIGHT_RIGHT"]
        # Save it on the properties list (one list of widths etc per road with all data from cells intersected by this street)
        street_index = streets_df[streets_df.STREET_ID == street_id].index[0]
        cells_list = street_properties_detailled['CELLS'][street_index]
        cell_index = cells_list.index(poly_num)
        street_properties_detailled['WIDTHS'][street_index][cell_index] = street_width
        street_properties_detailled['LEFT_HEIGHTS'][street_index][cell_index] = street_height_left
        street_properties_detailled['RIGHT_HEIGHTS'][street_index][cell_index] = street_height_right
    time_end = time.gmtime(time.time() - init)
    line = str(time_end.tm_hour) + 'hour' + str(time_end.tm_min) + 'min' + str(time_end.tm_sec)+ 'sec.'
    #print(poly_num, "poly  - time:", line)
time_end = time.gmtime(time.time() - init)
line = str(time_end.tm_mday) + "days"+ str(time_end.tm_hour) + 'hour' + str(time_end.tm_min) + 'min' + str(time_end.tm_sec)+ 'sec.'
print("All cells properties done  - time:", line)

# Average properties for streets intersecting multiple cells and update street width if the building method got it
urban_morphology_attributes = ['STREET_ID', 'WIDTH', 'LEFT_HEIGHT', 'RIGHT_HEIGHT', 'HEIGHT','LENGTH']
if save_canyon_attribute:
    urban_morphology_attributes.extend(['OPEN','RATIO'])
street_attributes = pd.DataFrame(columns=urban_morphology_attributes)
street_attributes['STREET_ID'] = streets_df['STREET_ID']
if save_canyon_attribute:
    street_attributes['OPEN'] = ['no']*n_streets
for street_index, is_tunnel in enumerate(streets_df['TUNNEL']) :
    # Height attribute
    if not is_tunnel:
        distribution = streets_df.loc[street_index,'DISTRIBUTION']
        left_heights = street_properties_detailled['LEFT_HEIGHTS'][street_index]
        right_heights = street_properties_detailled['RIGHT_HEIGHTS'][street_index]
        left_height = np.average(left_heights, weights=distribution)
        right_height = np.average(right_heights, weights=distribution)
        street_height = np.nanmean([left_height, right_height])
        street_attributes.loc[street_index,'LEFT_HEIGHT'] = left_height
        street_attributes.loc[street_index,'RIGHT_HEIGHT'] = right_height
        street_attributes.loc[street_index, 'HEIGHT'] = street_height
    else :
        street_attributes.loc[street_index, 'LEFT_HEIGHT'] = np.nan
        street_attributes.loc[street_index, 'RIGHT_HEIGHT'] = np.nan
        street_attributes.loc[street_index, 'HEIGHT'] = np.nan
    # Width attribute
    widths = street_properties_detailled['WIDTHS'][street_index]
    if np.isnan(np.min(widths)):  # if there is at least one Nan
        street_width = np.nanmean(widths)
    else : # if we found a width in every cell, compute the weighted average
        street_width = np.average(widths, weights=distribution)
    street_attributes.loc[street_index, 'WIDTH'] = street_width
    # Canyon street attribute (if the road is open on one or both sides)
    if save_canyon_attribute:
        if not is_tunnel:
            if (left_height > 0) and (right_height > 0) and (street_width > 0):
                street_attributes.loc[street_index, 'RATIO'] = street_height / street_width
            else:
                if (left_height == 0.0) and (right_height > 0) :
                    street_attributes.loc[street_index, 'OPEN'] = 'L'
                elif (right_height == 0.0) and (left_height > 0) :
                    street_attributes.loc[street_index, 'OPEN'] = 'R'
                elif (left_height == 0.0) and (right_height == 0.0) :
                    street_attributes.loc[street_index, 'OPEN'] = 'L/R'
                street_attributes.loc[street_index, 'RATIO'] = np.nan
        else :
            street_attributes.loc[street_index, 'RATIO'] = np.nan
    if street_index % 100 == 0:
        time_end = time.gmtime(time.time() - init)
        line = str(time_end.tm_hour) + 'hour' + str(time_end.tm_min) + 'min' + str(time_end.tm_sec)+ 'sec.'
        print(street_index, "streets  - time:", line)
time_end = time.gmtime(time.time() - init)
line = str(time_end.tm_mday) + "days"+ str(time_end.tm_hour) + 'hour' + str(time_end.tm_min) + 'min' + str(time_end.tm_sec)+ 'sec.'
print("All attributes synthetized. Now saving... Time:", line)

# Add canyon streets and distribution of streets in grid cells if necessary
attributes = ['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY','LENGTH','CELLS', 'DISTRIBUTION']
if not save_streets_on_grid :
    attributes.remove('CELLS')
    attributes.remove('DISTRIBUTION')
geometric_information = streets_df[attributes]
street_attributes = pd.merge(left=geometric_information, right=street_attributes, how='left', on=['STREET_ID'])

# Save road attributes
f = open(foutput + '.dat', 'wb')
pickle.dump(street_attributes, f)
f.close()
street_attributes.to_csv(foutput+'.csv',index=False)

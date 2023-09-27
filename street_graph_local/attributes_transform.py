##################################################################################
# Study: Street network generation
# Purpose: Attribute properties to streets in the Île-de-France graph
# Author: Marjolaine Lannes
# Creation date: February 28
# Note: Get height/width data from each cell and calculate each street weighted mean height/width
##################################################################################
import pickle
import pandas as pd
import numpy as np

# Parameters
aspect_ratio_lim = 1/5

# Input files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
f_streets = path + "temp/street_graph/graph_transform/street_links_width_distrib.dat"
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
street_height_f = path + 'temp/street_graph/attributes/streets_height_width.dat'
foutput = path + 'output/R_0_2/street_attributes.dat'
foutput_canyon = path + 'output/R_0_2/canyon_streets.dat'
foutput_csv = path + 'output/R_0_2/street_attributes.csv'

# Load data
streets_df = pickle.load(open(f_streets, 'rb'))
#streets_df = pd.DataFrame(streets_data)
# columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY', 'XA', 'XB', 'YA', 'YB', 'ROADS_WIDTH', 'SIDEWALKS_WIDTH', 'TOTAL_WIDTH', 'CELLS', 'DISTRIBUTION']
n_streets = streets_df.shape[0]
grid = pickle.load(open(gridfile,'rb'))
print('Data loaded...')

# Initialize data
street_attributes = pd.DataFrame(columns=['STREET_ID', 'WIDTH', 'HEIGHT','LEFT_HEIGHT', 'RIGHT_HEIGHT', 'CANYON','RATIO'])
street_attributes['STREET_ID'] = streets_df['STREET_ID']

# load data
street_heights_df = pickle.load(open(street_height_f, 'rb'))
# street_heights_df = pd.DataFrame(street_heights_data)
print(street_heights_df)
# Get street properties
for street_SUB_index, street_id in enumerate(list(street_heights_df['STREET_ID'])):
    # Street width/height
    street_width = street_heights_df.loc[street_SUB_index, "WIDTH"]
    street_height_left = street_heights_df.loc[street_SUB_index, "HEIGHT_LEFT"]
    street_height_right = street_heights_df.loc[street_SUB_index, "HEIGHT_RIGHT"]
    street_height = np.nanmean([street_height_left,street_height_right])
    # Save it on the properties list (one list of widths etc per road with all data from cells intersected by this street)
    street_index = streets_df[streets_df.STREET_ID == street_id].index[0]
    street_attributes.loc[street_SUB_index, 'HEIGHT'] = street_height
    street_attributes.loc[street_SUB_index,'LEFT_HEIGHT'] = street_height_left
    street_attributes.loc[street_SUB_index,'RIGHT_HEIGHT'] = street_height_right
    if np.isnan(street_width):
        matching_width = streets_df.loc[street_index,'TOTAL_WIDTH']
        street_attributes.loc[street_SUB_index,'WIDTH'] = matching_width
        street_width = matching_width
    else :
        street_attributes.loc[street_SUB_index,'WIDTH'] = street_width
        #if street_width > 0 :
            #print("Calculated width ok")
        #else :
            #print("No calculated width")
    if (street_height_left > 0) and (street_height_right > 0) and (street_width > 0) :
        coeff = street_height / street_width
        street_attributes.loc[street_SUB_index, 'RATIO'] = coeff
        if coeff > aspect_ratio_lim :
            street_attributes.loc[street_SUB_index, 'CANYON'] = True
        else :
            print("could have been", coeff, street_height, street_width)
            street_attributes.loc[street_SUB_index, 'CANYON'] = False
    else :
        street_attributes.loc[street_SUB_index, 'CANYON'] = False
        street_attributes.loc[street_SUB_index, 'RATIO'] = np.nan

geometric_information = streets_df[['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY']]
street_attributes = pd.merge(left=geometric_information, right=street_attributes, how='left', on=['STREET_ID'])
canyon_streets = street_attributes[street_attributes.CANYON == True]
print(canyon_streets)

# Save road attributes
pickle.dump(street_attributes, open(foutput, 'wb'))
street_attributes.to_csv(foutput_csv, index=False)
pickle.dump(canyon_streets, open(foutput_canyon, 'wb'))

##################################################################################
# Study: Street network generation
# Purpose: Get streets attributes: width, height, cells distribution
# Author: Marjolaine Lannes
# Creation date: April 25, 2023
# Note: # penser au traitement des cross roads du point de vue des émissions (pas pris en compte dans le calcul de nodes)
##################################################################################
import pickle
import pandas as pd
import geopandas as gpd
from itertools import repeat
import time

# Input and output files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
sidewalks_f = path + "temp/street_graph/attributes/Paris_sidewalks.dat"
# buildings_f = path + 'temp/load_BDTOPO/buildings.dat'
nodes_f = path + 'temp/street_graph/groups/nodes.dat'
roads_widths_dir = path + "temp/road_graph/subdomains_roads_widths/"
transform_dir = path + "temp/street_graph/graph_transform/"
street_nodes_f = transform_dir + "street_nodes.csv"
street_links_f = transform_dir + "street_links.dat"
links_transform_f = transform_dir + "links_transform.dat"
nodes_transform_f = transform_dir + "nodes_transform.dat"
links_inverse_transform_f = transform_dir + "links_inverse_transform.dat"
street_links_attributes_f = transform_dir + "street_links_width_distrib.dat"
street_cells_distrib_f = transform_dir + "street_cells_distrib.csv"
grid_file = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
outdir = transform_dir + 'street_cells/'

# Get cell information
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
grid=pickle.load(open(gridfile,'rb'))
poly_num = 4672 #int(sys.argv[1])
x0 = grid['we_id'][poly_num]
y0 = grid['sn_id'][poly_num]
prefix = 'we_id_from_' + str(x0) + '_sn_id_from_' + str(y0)
buildings_dir = path + 'temp/load_BDTOPO/subdomains_buildings_with_buffer_100m_lamb93/'
fname = 'sub_' + prefix + '.dat'
fname2 = 'sub_' + 'we_id_from0_' + str(x0) + '_sn_id_from0_' + str(y0) + '.dat'

# Load data
Street_links = pickle.load(open(street_links_f,"rb")) # columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY', 'XA', 'XB', 'YA', 'YB']
Street_nodes = pd.read_csv(street_nodes_f) # columns = ['NODE_ID', 'X', 'Y']
links_transform = pickle.load(open(links_transform_f,"rb"))
# columns=['ROAD_ID', 'GROUP_ID', 'PARALLEL', 'MERGED', 'COLOR', 'GEOMETRY','STREET_IDS', 'SHARES']
nodes_transform = pickle.load(open(nodes_transform_f,"rb")) # columns=['ROAD_NODE', 'GROUP_NODE', 'STREET_NODE']
sw_f = open(sidewalks_f, 'rb')
sidewalks = pickle.load(sw_f)
sw_f.close()
buildings = pickle.load(open(buildings_dir + fname, 'rb'))
# buildings = pickle.load(open(buildings_f, 'rb'))
roads_widths = pickle.load(open(roads_widths_dir + fname2, 'rb'))

# 1) Inverse transform : streets to roads ids
links_inverse_transform = {"STREET_ID":[],"ROADS_IDS":[],'ROADS_SHARES':[]}
streets_list = list(Street_links["STREET_ID"])
links_inverse_transform["STREET_ID"] = streets_list
n_streets = len(streets_list)
links_inverse_transform['ROADS_IDS'] = [[] for i in repeat(None, n_streets)]
links_inverse_transform['ROADS_SHARES'] = [[] for j in repeat(None, n_streets)]
for r, road_r in enumerate(links_transform["ROAD_ID"]):
    matching_streets = links_transform["STREET_IDS"][r]
    shares = links_transform["SHARES"][r]
    for s, street_s in enumerate(matching_streets) :
        street_index = streets_list.index(street_s)
        links_inverse_transform["ROADS_IDS"][street_index].append(road_r)
        links_inverse_transform["ROADS_SHARES"][street_index].append(shares[s])
print("Inverse link transform computed...")

# 2) Widths ## IF NO BUILDING FOUND AT LEFT / RIGHT SIDE
roads_list = links_transform['ROAD_ID']
for i, street_id in enumerate(streets_list):
    # for each street, sum parallel roads widths
    matching_roads = links_inverse_transform['ROADS_IDS'][i]
    parallel_roads_widths = []
    for road_id in matching_roads :
        road_index = roads_list.index(road_id)
        if links_transform['PARALLEL'][road_index] == True :
            parallel_roads_widths.append(roads_widths['mean_width'][road_index])
    width = sum(parallel_roads_widths)
    # add sidewalks widths :
    n_sidewalks = len(sidewalks['SIDEWALKS'][i])
    if n_sidewalks > 0 : # calculate sw_width
        sw_width = n_sidewalks * float(sidewalks['MEAN_SW_WIDTH'][i])
    else : # If no sidewalk is registered, calculate sw_width based on the street width
        n_sidewalks = 2
        if width <= 4.0 :
            sw_width = n_sidewalks * 1.77 # Q1 sidewalk width
        elif 4.0 < width <= 5.0 : # Q1 and Q3 of width
            sw_width = n_sidewalks * 4.37 # mean sidewalk width
        else :
            sw_width = n_sidewalks * 5.83 # Q3 sidewalk width
    # save total street width
    Street_links.loc[i, 'ROADS_WIDTH'] = width
    Street_links.loc[i, 'SIDEWALKS_WIDTH'] = sw_width
    Street_links.loc[i, 'TOTAL_WIDTH'] = width + sw_width
print("Street widths computed...")

# Save data
links_inverse_transform_df = pd.DataFrame.from_dict(links_inverse_transform)
links_inverse_transform_df.to_csv(links_inverse_transform_f)
# Street_links.to_csv(street_links_attributes_f)
pickle.dump(Street_links, open(street_links_attributes_f, "wb"))

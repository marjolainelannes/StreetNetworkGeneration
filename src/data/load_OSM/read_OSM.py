##################################################################################
# Study: Street network generation
# Purpose: Create a graph for the road network
# Note: Read edges/nodes attributes of OSM road network
##################################################################################
import pandas as pd
import xml.etree.ElementTree as Xet
import os

# Config
path = "../../../"
input_dir = path + "data/"
cache_dir  = path + "temp/"

# files
MATSim_network = input_dir + "OSM/ile_de_france_network.xml"
output_links_f = cache_dir + "load_OSM/OSM_links.csv"
output_nodes_f = cache_dir + "load_OSM/OSM_nodes.csv"
if not os.path.exists(cache_dir + "load_OSM/"):
    os.makedirs(cache_dir + "load_OSM/")

# Load OSM information
xmlparse = Xet.parse(MATSim_network)
root = xmlparse.getroot()
nodes = root[0]
n_nodes = len(nodes)
links = root[1]
n_links = len(links)

# Generate road_nodes dataframe
nodes_data = []
for i in range(0,n_nodes):
    node_i = nodes[i]
    row_i = node_i.attrib
    nodes_data.append(row_i)
nodes_df = pd.DataFrame(nodes_data)

# Generate road_links dataframe
links_data = []
for i in range (0, n_links) : #one dataframe for all elements of this type
    road_i = links[i]
    # Get the road default attributes: from, to, length, etc
    row_i = road_i.attrib
    # Add 'roundabout' and 'tunnel' attributes
    row_i['roundabout'] = False
    row_i['tunnel'] = False
    for attributes in road_i.iter('attributes'):
         for attribute in attributes.iter('attribute') :
             attribute_name = attribute.get('name')
             # check if the road is a roundabout
             if attribute_name == "osm:way:junction" :
                 junction_type = attribute.text
                 if junction_type in ['roundabout', 'circular'] : # it is also possible to add jughandle
                     row_i['roundabout'] = True
             # check if the road is a tunnel
             if attribute_name == "osm:way:tunnel" :
                 row_i['tunnel'] = True
    links_data.append(row_i)
links_df = pd.DataFrame(links_data)

# Save data
links_df.to_csv(output_links_f, index = False)
nodes_df.to_csv(output_nodes_f, index = False)

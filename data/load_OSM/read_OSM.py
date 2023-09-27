##################################################################################
# Study: Street network generation
# Purpose: Create a geograph for the road network
# Author: Marjolaine Lannes
# Creation date: May 5, 2023
# Note: Read edges/nodes attributes of OSM road network
##################################################################################
import pandas as pd
import xml.etree.ElementTree as Xet

# Input files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
MATSim_network = path + "data/OSM/ile_de_france_network.xml"
output_links_f = path + "temp/load_OSM/OSM_links.csv"
output_nodes_f = path + "temp/load_OSM/OSM_nodes.csv"

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
    # Add attributes 'roundabout', 'Boulevard_peripherique' and 'tunnel'
    row_i['roundabout'] = False
    row_i['tunnel'] = False
    for attributes in road_i.iter('attributes'):
         for attribute in attributes.iter('attribute') :
             attribute_name = attribute.get('name')
             # check if the road is a roundabout
             if attribute_name == "osm:way:junction" :
                 junction_type = attribute.text
                 if junction_type in ['roundabout', 'circular'] : # also : jughandle
                     row_i['roundabout'] = True
             # check if the road is a tunnel
             if attribute_name == "osm:way:tunnel" :
                 row_i['tunnel'] = True
    links_data.append(row_i)
links_df = pd.DataFrame(links_data)

# Save data
links_df.to_csv(output_links_f, index = False)
nodes_df.to_csv(output_nodes_f, index = False)

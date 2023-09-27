# To transform emissions on roads in emissions on streets :
# calculate the share of a road's emissions going to a street

import pickle
import pandas as pd

# Input files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
links_transform_f = path + "temp/street_graph/graph_transform/" + "links_transform.dat"
roads_preproc_f = path + "temp/load_OSM/roads_preprocessing.csv"
f_canyon = path + 'output/R_0_2/canyon_streets.dat'
output_f = path + "output/R_0_2/roads_to_streets_canyon.csv"

# Load data
links_transform = pickle.load(open(links_transform_f,'rb'))
# dict : ROAD_ID, STREETS_IDS, SHARES ( + groups)
roads_preproc = pd.read_csv(roads_preproc_f)
# df : ['link_id','matching_link','pt']
canyons = pickle.load(open(f_canyon, "rb"))
canyons_list = list(canyons['STREET_ID'])
print(canyons_list)
Roads_to_Street = pd.DataFrame(columns=['ROAD_ID', 'STREET_ID','SHARE'])
roads_after_filter_list = list(links_transform['ROAD_ID'])
print(links_transform.keys())

# Get the transform
for i, road_id in enumerate(roads_after_filter_list) :
    matching_streets = links_transform['STREET_IDS'][i]
    n_streets = len(matching_streets)
    for j in range(n_streets):
        street = matching_streets[j]
        if street in canyons_list :
            new_row = pd.DataFrame(columns=['ROAD_ID', 'STREET_ID','SHARE'])
            new_row.loc[0,'ROAD_ID'] = road_id
            new_row.loc[0, 'STREET_ID'] = street
            new_row.loc[0, 'SHARE'] = links_transform['SHARES'][i][j]
            Roads_to_Street = pd.concat([Roads_to_Street,new_row], ignore_index=True)

# Get the transform for roads that we merged during preproc
for i, road_id in enumerate(roads_preproc['link_id']) :
    updated_road = roads_preproc.loc[i, 'matching_link']
    if (updated_road in roads_after_filter_list) and (str(roads_preproc.loc[i,'pt']) == "False") :
        updated_road_index = roads_after_filter_list.index(updated_road)
        matching_streets = links_transform['STREET_IDS'][updated_road_index]
        n_streets = len(matching_streets)
        for j in range(n_streets):
            street = matching_streets[j]
            if street in canyons_list:
                new_row = pd.DataFrame(columns=['ROAD_ID', 'STREET_ID', 'SHARE'])
                new_row.loc[0, 'ROAD_ID'] = road_id
                new_row.loc[0, 'STREET_ID'] = street
                new_row.loc[0, 'SHARE'] = links_transform['SHARES'][updated_road_index][j]
                Roads_to_Street = pd.concat([Roads_to_Street, new_row], ignore_index=True)

# Save data
Roads_to_Street.to_csv(output_f, index=False)

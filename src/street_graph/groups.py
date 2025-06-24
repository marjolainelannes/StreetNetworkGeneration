##################################################################################
# Study: Street network generation
# Purpose: To identify roads to merge within a street
##################################################################################
import os, time, pickle, sys
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Polygon,Point
import numpy as np
sys.path.append("../../include")
from lines_geometry import *
from create_poly_around_road import create_1_poly_around_road

# Parameters: thresholds for merging
n_cells = 8250
links_distance = 500 # meters
poly_width = 100  # meters
building_area_threshold = 0.05
road_length_percentage = 0.8
number_of_colors = 25
building_height_threshold = 3.0 # meters, 10% of buildings in the region
parallels_angle = 30 # degrees

# Directories
path = "../../"
data_dir = path + "data/"
cache_dir = path + "temp/"
osm_dir = cache_dir + "data/load_OSM/"
groups_dir = cache_dir + 'street_graph/groups/'
if not os.path.exists(groups_dir):
    os.makedirs(groups_dir)
  
# Input and output files
adjacent_cells_f = data_dir + "grid_data/adjacent_cells.dat"
grid_file = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
roads_f = osm_dir + "roads_links.csv"
nodes_f = osm_dir + "roads_nodes.csv"
detailed_roads_f = osm_dir + "OSM_links.csv"
neighborhood_links_f = osm_dir + "neighborhood_links_" + str(links_distance)+ "m.dat"
roads_per_cell_f = osm_dir + "roads_list_per_cell.dat"
buildings_dir = path + 'temp/data/urban_form/subdomains_buildings_with_adjacent_cells_lamb93/'
waterways_f = osm_dir + "waterways/filtered_waterways.shp"
fout_nodes = groups_dir + 'group_nodes.dat'
fout_roads = groups_dir + 'roads_to_groups.csv'
fout_roads_dat = groups_dir + 'roads_to_groups.dat'
fout_groups = groups_dir + 'groups.dat'
init_f = groups_dir + 'roads_to_groups_init.csv'

# Load and format data
adjacent_cells = pickle.load(open(adjacent_cells_f, "rb"))
grid = pickle.load(open(grid_file, 'rb'))
roads = pd.read_csv(roads_f)
roads_per_cell = pickle.load(open(roads_per_cell_f, "rb"))
roads_ids = [str(road) for road in roads['ID']]
roads['ROAD_ID'] = roads_ids
roads['GEOMETRY'] = roads['GEOMETRY'].apply(lambda x: shapely.wkt.loads(x))
detailed_roads = pd.read_csv(detailed_roads_f)
detailed_roads = detailed_roads.astype({'id':'str'})
neighborhood_links = pickle.load(open(neighborhood_links_f, 'rb'))
nodes_types = {'ID':str,'X':float, 'Y':float, 'LINKS':object}
nodes_OSM = pd.read_csv(nodes_f, dtype=nodes_types)
nodes_OSM['LINKS'] = nodes_OSM['LINKS'].apply(eval)
waterways = gpd.read_file(waterways_f)
list_of_cells = np.arange(0,n_cells)
f = open(groups_dir+'groups_run.txt', 'w')
print('Data loaded... ', time.time())

# Initialize data
Roads_to_groups = pd.DataFrame(columns=['ROAD_ID', 'GROUP_ID', 'MAIN_ROAD', 'MERGED', 'COLOR', 'ROUNDABOUT', 'TUNNEL', 'GEOMETRY'])
detailed_roads_list = list(detailed_roads['id'])
for i, road_i in enumerate(roads_ids) :
    Roads_to_groups.loc[i, 'ROAD_ID'] = road_i
    # Get the road detailed information (if it is a roundabout, a tunnel, its geometry)
    road_index = detailed_roads_list.index(road_i)
    Roads_to_groups.loc[i,'ROUNDABOUT'] = detailed_roads.loc[road_index,'roundabout']
    Roads_to_groups.loc[i,'TUNNEL'] = detailed_roads.loc[road_index,'tunnel']
    Roads_to_groups.loc[i,'GEOMETRY'] = roads.loc[i,'GEOMETRY']
Roads_to_groups = Roads_to_groups.astype({'ROAD_ID':str})
Roads_to_groups.to_csv(init_f, index=False)
print("Roads_to_groups", Roads_to_groups.shape[0])
print("roads", roads.shape[0])

# Filter roundabouts
Roundabouts = Roads_to_groups[Roads_to_groups.ROUNDABOUT == True]
roundabouts_ids = list(Roundabouts['ROAD_ID'])

# 1. Identify parallel roads
## purpose: to test which roads have to be added to the group of roads to merge with road_i, if there is
print('1. Identify parallel roads')
f.write('1. Identify parallel roads\n')
# 1.a) First identification
start_parallel = time.time()
roads_of_parallels_id = {}
parallels_nodes = {}
longest_road_within_parallels = {}
parallels_max_length = {}
tunnel = {}
parallels_per_cell = {}
parallel_id_of_road = {}
non_merged_roads = [item for item in roads_ids if (item not in roundabouts_ids)]
n_parallels = 0
for cell_id in list_of_cells :
    start_cell = time.time()
    roads_cell_i = roads_per_cell[cell_id]
    roads_ids_cell_without_roundabouts = [item for item in roads_cell_i if (item in non_merged_roads)]
    if len(roads_ids_cell_without_roundabouts) == 0 :
        parallels_per_cell[cell_id] = []
        line = "Step 1 parallel roads - Cell" + str(cell_id) + ", number of parallel roads :" + str(len(longest_road_within_parallels))
        print(line)
        f.write(line + '\n')
        continue
    x = grid['we_id'][cell_id]
    y = grid['sn_id'][cell_id]
    buildings_f = buildings_dir + 'sub_we_id_from0_' + str(x) + '_sn_id_from0_' + str(y) + '.dat'
    buildings_data = pickle.load(open(buildings_f, 'rb'))
    buildings_df = pd.DataFrame(buildings_data)
    buildings = gpd.GeoDataFrame(data=buildings_df, geometry=buildings_df['POLYGON'])
    buildings = buildings[buildings.HEIGHT > building_height_threshold]
    parallels_cell_i = []
    for road_i in roads_ids_cell_without_roundabouts :
        i = roads_ids.index(road_i)
        if road_i in non_merged_roads :
            # Get road geometry
            node_A = roads.loc[i,'NODE_A']
            node_B = roads.loc[i,'NODE_B']
            road_xa = roads.loc[i,'XA']
            road_ya = roads.loc[i,'YA']
            road_xb = roads.loc[i,'XB']
            road_yb = roads.loc[i,'YB']
            node_A_point = Point(road_xa, road_ya)
            points_coordinates_1 = [road_xa, road_ya, road_xb, road_yb]
            length_i = roads.loc[i,"LENGTH"]
            tunnel_i = Roads_to_groups.loc[i, 'TUNNEL']
            poly_i = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, poly_width)
            # Initialize merging
            non_merged_roads.remove(road_i)
            merged = False
            n_parallels = len(longest_road_within_parallels)
        else :
            merged = True
        roads_representing_parallels = list(longest_road_within_parallels.values()) + non_merged_roads  # add all roads that are not merged for the moment
        roads_to_test = neighborhood_links[road_i] # test only links in neighborhood
        cnt = 0 # counter of the place within roads_to_test
        while not (merged or cnt==len(roads_to_test)):
            road_j = roads_to_test[cnt]
            j = roads_ids.index(road_j)
            # filter if tunnel or not
            tunnel_j = Roads_to_groups.loc[j, 'TUNNEL']
            if (tunnel_i != tunnel_j) or (road_j in roundabouts_ids) or (i==j):
                cnt = cnt + 1
                continue
            # data for tests
            road_xc = roads.loc[j,'XA']
            road_yc = roads.loc[j,'YA']
            road_xd = roads.loc[j,'XB']
            road_yd = roads.loc[j,'YB']
            node_C_point = Point(road_xc, road_yc)
            node_D_point = Point(road_xd, road_yd)
            node_C = roads.loc[j,'NODE_A']
            node_D = roads.loc[j,'NODE_B']
            edge_j = roads.loc[j,"GEOMETRY"]
            length_j = roads.loc[j,"LENGTH"]
            points_coordinates_2 = [road_xc, road_yc, road_xd, road_yd]
            if poly_i.intersects(edge_j):
                intersect_i_j = poly_i.intersection(edge_j).length
            else:
                intersect_i_j = 0
            # test road_i vs road_j
            if ((intersect_i_j > road_length_percentage * length_i) or (intersect_i_j > road_length_percentage * length_j)) and parallel_lines(points_coordinates_1, points_coordinates_2, parallels_angle) :
                # Check that there is no building between the two roads (condition on area)
                if shapely.distance(node_A_point, node_C_point) > shapely.distance(node_A_point, node_D_point) :
                    poly_roads = Polygon(((road_xa, road_ya),(road_xb, road_yb),(road_xc, road_yc),(road_xd, road_yd)))
                else :
                    poly_roads = Polygon(((road_xa, road_ya), (road_xb, road_yb), (road_xd, road_yd), (road_xc, road_yc)))
                subdata = buildings[buildings.geometry.intersects(poly_roads)]
                #subdata = subdata[subdata.HEIGHT > building_height_threshold]
                total_buildings_area = sum(subdata.geometry.area)
                intersecting_waterways = waterways[waterways.geometry.intersects(poly_roads)]
                if (total_buildings_area < building_area_threshold * poly_roads.area) and (
                        sum(intersecting_waterways.geometry.length) == 0):   # if there is no building large enough, then merge the two roads
                    if road_j in roads_representing_parallels:
                        road_j_rep = road_j
                    else:
                        parallels_j = parallel_id_of_road[road_j]
                        road_j_rep = longest_road_within_parallels[parallels_j]
                    parallels_cnt = roads_representing_parallels.index(road_j_rep)
                    if (n_parallels > 0) and (parallels_cnt + 1 <= n_parallels) :  # if road_j is already in a set of parallels, add road_i to it
                        parallels_id = list(longest_road_within_parallels.keys())[parallels_cnt]
                        roads_of_parallels_id[parallels_id].append(road_i)
                        parallels_nodes[parallels_id].extend([node_A, node_B])
                        parallel_id_of_road[road_i] = parallels_id
                        if length_i > length_j:
                            parallels_max_length[parallels_id] = length_i
                            longest_road_within_parallels[parallels_id] = road_i
                    else:  # else, create a new set of parallels with road_i and road_j
                        parallels_id = n_parallels
                        n_parallels += 1
                        parallel_id_of_road[road_i] = parallels_id
                        parallel_id_of_road[road_j] = parallels_id
                        roads_of_parallels_id[parallels_id] = [road_i, road_j]
                        parallels_nodes[parallels_id] = [node_A, node_B, node_C, node_D]
                        parallels_cell_i.append(parallels_id)
                        tunnel[parallels_id] = tunnel_i
                        if road_j in non_merged_roads:
                            non_merged_roads.remove(road_j)
                        if length_j > length_i:
                            longest_road_within_parallels[parallels_id] = road_j
                            parallels_max_length[parallels_id] = length_j
                        else:
                            longest_road_within_parallels[parallels_id] = road_i
                            parallels_max_length[parallels_id] = length_i
                    merged = True
                else :
                    cnt = cnt + 1
            else :
                cnt = cnt + 1
        if not merged :
            non_merged_roads.append(road_i)
    parallels_per_cell[cell_id] = parallels_cell_i
    line = "Step 1 parallel roads - Cell" + str(cell_id) +  ", number of parallel roads :"+ str(len(longest_road_within_parallels))
    print(line)
    f.write(line+'\n')
    time_cell_sec = time.time() - start_cell
    time_cell = time.gmtime(time_cell_sec)
    line = str(time_cell.tm_mday )+ "days"+ str(time_cell.tm_hour) + 'hour' + str(time_cell.tm_min) + 'min' + str(time_cell.tm_sec)+ 'sec.'
    print(line)
    f.write(line + '\n')
parallels_ids = list(np.arange(0,n_parallels))
time_parallell_roads_sec = time.time()-start_parallel
time_PR = time.gmtime(time_parallell_roads_sec)
line = "1.a. parallels, took" + str(time_PR.tm_mday )+ "days" + str(time_PR.tm_hour) + 'hour' + str(time_PR.tm_min) + 'min' + str(time_PR.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')
# In the end : non_merged_roads and parallels with groups of roads

# 1.b) Check if 2 pairs of parallel roads are parallel to each other
line= "1.b) Check if 2 pairs of parallel roads are parallel to each other"
print(line)
f.write(line + '\n')
## then the 2 sets are to be merged, based on longest road of each group
roads_to_test = list(longest_road_within_parallels.values())
for cell_id in list_of_cells :
    print(cell_id, "cell: merge parallels")
    x = grid['we_id'][cell_id]
    y = grid['sn_id'][cell_id]
    buildings_f = buildings_dir + 'sub_we_id_from0_' + str(x) + '_sn_id_from0_' + str(y) + '.dat'
    buildings_data = pickle.load(open(buildings_f, 'rb'))
    buildings_df = pd.DataFrame(buildings_data)
    buildings = gpd.GeoDataFrame(data=buildings_df, geometry=buildings_df['POLYGON'])
    buildings = buildings[buildings.HEIGHT > building_height_threshold]
    parallels_to_test = []
    for cell_j in adjacent_cells[cell_id]:
        if cell_j in parallels_per_cell.keys():
            parallels_to_test = list(set(parallels_to_test + parallels_per_cell[cell_j]))
    final_parallels_cell_i = parallels_per_cell[cell_id]
    for parallels_i in parallels_per_cell[cell_id] :
        road_i = roads_of_parallels_id[parallels_i][0]
        i = roads_ids.index(road_i)
        # purpose: to test which roads have to be added to the set of roads to merge with road_i, if there is
        if road_i in roads_to_test :
            # Get road geometry
            node_A = roads.loc[i, 'NODE_A']
            node_B = roads.loc[i, 'NODE_B']
            road_xa = roads.loc[i, 'XA']
            road_ya = roads.loc[i, 'YA']
            road_xb = roads.loc[i, 'XB']
            road_yb = roads.loc[i, 'YB']
            node_A_point = Point(road_xa, road_ya)
            points_coordinates_1 = [road_xa, road_ya, road_xb, road_yb]
            length_i = roads.loc[i, "LENGTH"]
            tunnel_i = tunnel[parallels_i]
            # Define polygon
            poly_i = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, poly_width)
            # Initialize merging
            merged = False
        else :
            merged = True
        cnt = 0 # counter of the place within roads_to_test
        while not (merged or cnt == len(parallels_to_test)):
            parallels_j = parallels_to_test[cnt]
            if (parallels_j == parallels_i) or (parallels_j not in parallels_ids) :
                cnt += 1
                continue
            # filter tunnel (or not)
            tunnel_j = tunnel[parallels_j]
            if tunnel_i != tunnel_j :
                cnt += 1
                continue
            road_j = longest_road_within_parallels[parallels_j] 
            j = roads_ids.index(road_j)
            # data for tests
            road_xc = roads.loc[j,'XA']
            road_yc = roads.loc[j,'YA']
            road_xd = roads.loc[j,'XB']
            road_yd = roads.loc[j,'YB']
            node_C_point = Point(road_xc, road_yc)
            node_D_point = Point(road_xd, road_yd)
            node_C = roads.loc[j,'NODE_A']
            node_D = roads.loc[j,'NODE_B']
            edge_j = roads.loc[j,"GEOMETRY"]
            length_j = roads.loc[j,"LENGTH"]
            points_coordinates_2 = [road_xc, road_yc, road_xd, road_yd]
            if poly_i.intersects(edge_j):
                intersect_i_j = poly_i.intersection(edge_j).length
            else:
                intersect_i_j = 0
            ## Test: poly_around_road
            if ((intersect_i_j > road_length_percentage * length_i) or (intersect_i_j > road_length_percentage * length_j)) and parallel_lines(points_coordinates_1, points_coordinates_2, parallels_angle): # then we only keep group_j and merge it with group_i
                # Check that there is no building between the two roads (condition on area)
                if shapely.distance(node_A_point, node_C_point) > shapely.distance(node_A_point, node_D_point) :
                    poly_roads = Polygon(((road_xa, road_ya), (road_xb, road_yb), (road_xc, road_yc), (road_xd, road_yd)))
                else:
                    poly_roads = Polygon(((road_xa, road_ya), (road_xb, road_yb), (road_xd, road_yd), (road_xc, road_yc)))
                subdata = buildings[buildings.geometry.intersects(poly_roads)]
                # subdata = subdata[subdata.HEIGHT > building_height_threshold]
                total_buildings_area = sum(subdata.geometry.area)
                intersecting_waterways = waterways[waterways.geometry.intersects(poly_roads)]
                if (total_buildings_area < building_area_threshold * poly_roads.area) and (
                        sum(intersecting_waterways.geometry.length) == 0): # if there is no building large enough, then merge the two roads
                    # adding road_i in parallel_roads to the corresponding list
                    if length_i > length_j:
                        longest_road_within_parallels[parallels_j] = road_i
                        parallels_max_length[parallels_j] = length_i
                        if road_j in roads_to_test :
                            roads_to_test.remove(road_j)
                    else:
                        if road_i in roads_to_test:
                            roads_to_test.remove(road_i)
                    merged = True
                    roads_of_parallels_id[parallels_j] = list(set(roads_of_parallels_id[parallels_i] + roads_of_parallels_id[parallels_j]))
                    parallels_nodes[parallels_j] = list(set(parallels_nodes[parallels_i] + parallels_nodes[parallels_j]))
                    roads_of_parallels_id.pop(parallels_i)
                    parallels_nodes.pop(parallels_i)
                    final_parallels_cell_i.remove(parallels_i)
                    parallels_ids.remove(parallels_i)
                    parallels_max_length.pop(parallels_i)
                    longest_road_within_parallels.pop(parallels_i)
                    tunnel.pop(parallels_i)
                else :
                    cnt += 1
            else:
                cnt += 1
    parallels_per_cell[cell_id] = final_parallels_cell_i
time_total = time.time() - start_parallel
time_t = time.gmtime(time_total)
line = "It took" + str(time_t.tm_mday )+ "days"+ str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

line = "2. Identify junctions among parallel roads"
print(line)
f.write(line + '\n')
start_2 = time.time()
parallels_junctions = pd.DataFrame(columns = ['PARALLELS_ID', 'JUNCTION_TO_PARALLELS'])
for cell_id in list_of_cells : # list_of_cells
    print(cell_id, "cell: roads = junction and parallel")
    parallels_to_test = []
    for cell_j in adjacent_cells[cell_id]:
        if cell_j in parallels_per_cell.keys():
            parallels_to_test = list(set(parallels_to_test + parallels_per_cell[cell_j]))
    ## First, for each road of a set of parallels i, check if it is a junction for roads in another set of parallels j (among parallels_to_test)
    for parallels_i in parallels_per_cell[cell_id] :
        roads_within_subset = roads_of_parallels_id[parallels_i]
        parallels_to_check = parallels_ids.copy()
        for road_i in roads_within_subset :
            i = roads_ids.index(road_i)
            node_A = str(roads.loc[i,'NODE_A'])
            node_B = str(roads.loc[i,'NODE_B'])
            cnt = 0
            # for each set of parallels to merge, check if road_i intersects at least two of the parallels' roads
            while cnt < len(parallels_to_test):
                parallels_j = parallels_to_test[cnt]
                # test junction
                if (parallels_j in parallels_to_check) and (parallels_i != parallels_j) and (tunnel[parallels_i] == tunnel[parallels_j]) :
                    nodes = [node_A, node_B]  # list of nodes that have not been identified as connected to the parallels yet
                    for node in parallels_nodes[parallels_j]:
                        if (str(node) == node_A) and (node_A in nodes):
                            nodes.remove(node_A)
                        elif (str(node) == node_B) and (node_B in nodes):
                            nodes.remove(node_B)
                    if len(nodes) == 0:  # i.e. if both nodes are connected to roads of the subset
                        parallels_i_j = [[parallels_i, parallels_j]]
                        row = pd.DataFrame(data = parallels_i_j, columns=['PARALLELS_ID', 'JUNCTION_TO_PARALLELS'])
                        parallels_junctions = pd.concat([parallels_junctions, row], ignore_index=True)
                        parallels_to_check.remove(parallels_j)
                    else:
                        cnt += 1
                else :
                    cnt += 1
## Then, save this information in two dictionnaries
parallels_with_junctions = {}
junction_to_parallels = {}
for p in parallels_ids :
    roads_in_p_are_junctions = parallels_junctions[parallels_junctions.PARALLELS_ID == p] # select lines where p has a // road which is a junction
    if roads_in_p_are_junctions.shape[0] > 0 :
        junction_to_parallels[p] = list(roads_in_p_are_junctions['JUNCTION_TO_PARALLELS']) # return the list of groups in which it is a junction
    parallels_with_junctions_of_p = parallels_junctions[parallels_junctions.JUNCTION_TO_PARALLELS == p]
    if parallels_with_junctions_of_p.shape[0] > 0 :
        parallels_with_junctions[p] = list(parallels_with_junctions_of_p['PARALLELS_ID']) # return the list of groups in which it is a junction

line = "2.b) Identify parallels such as p1 has a junction in p2 and vice-versa"
print(line)
f.write(line + '\n')
junction_to_parallels_init = junction_to_parallels.copy()
for p1, p2_list in junction_to_parallels_init.items() :
    for p2 in p2_list :
        if p2 in junction_to_parallels.keys() and p1 in junction_to_parallels[p2] :
            # check which one is the longest
            if parallels_max_length[p1] < parallels_max_length[p2] : # if it is p2,
                junction_to_parallels[p2].remove(p1) # then p1 has no more junction from p2
                parallels_with_junctions[p1].remove(p2)
                if len(junction_to_parallels[p2]) == 0:
                    junction_to_parallels.pop(p2)
                if len(parallels_with_junctions[p1]) == 0:
                    parallels_with_junctions.pop(p1)
            else : # vice-versa
                junction_to_parallels[p1].remove(p2)
                parallels_with_junctions[p2].remove(p1)
                if len(junction_to_parallels[p1]) == 0:
                    junction_to_parallels.pop(p1)
                if len(parallels_with_junctions[p2]) == 0:
                    parallels_with_junctions.pop(p2)
time_total = time.time() - start_2
time_t = time.gmtime(time_total)
line = "It took" + str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

line = "3. Deduce the structure of groups"
print(line)
f.write(line + '\n')
start_3 = time.time()
line = "3.a) Associate parallels within the subset of main roads or the subset of junctions of a group"
print(line)
f.write(line + '\n')
n_parallels = len(parallels_ids)
order = []
for p in parallels_ids : ## first, order parallels
    if p in junction_to_parallels.keys():
        order.append(1)
    else :
        order.append(0)
ordered_parallels_ids = [i for _, i in sorted(zip(order, parallels_ids))]
associated_groups = {"Main_roads":{}, "Junctions":{}}
n_groups=0
identified_parallels = []
cnt = 0
while cnt < n_parallels : # is it already in associated_groups ?
    parallels_i = ordered_parallels_ids[cnt]
    cnt += 1
    if parallels_i in identified_parallels :
        continue
    if parallels_i in junction_to_parallels.keys(): # if p_i has crossroads in p_j
        list_of_parallels_j = junction_to_parallels[parallels_i]
        existing_group = False
        j=0
        while not (existing_group or j==len(list_of_parallels_j)):
            parallels_j = list_of_parallels_j[j]
            j += 1
            if parallels_j in identified_parallels : # if p_j already identified in a group
                found = False
                cnt_group_j = 0
                while not found :
                    if cnt_group_j in associated_groups["Main_roads"].keys():
                        if parallels_j in associated_groups["Main_roads"][cnt_group_j] :
                            first_identified_parallels = associated_groups["Junctions"][cnt_group_j]
                            first_identified_parallels.append(parallels_i)
                            associated_groups["Junctions"][cnt_group_j] = first_identified_parallels
                            found = True
                        elif parallels_j in associated_groups["Junctions"][cnt_group_j] :
                            first_identified_parallels = associated_groups["Main_roads"][cnt_group_j]
                            first_identified_parallels.append(parallels_i)
                            associated_groups["Main_roads"][cnt_group_j] = first_identified_parallels
                            found = True
                        if found :
                            identified_parallels.append(parallels_i)
                            existing_group = True
                        else:
                            cnt_group_j += 1
                    else :
                        cnt_group_j += 1
        if not existing_group :
            parallels_j = list_of_parallels_j[0]
            associated_groups["Main_roads"][n_groups] = [parallels_j]
            junctions_subset_new_group = [parallels_k for parallels_k in parallels_with_junctions[parallels_j] if parallels_k not in identified_parallels]
            associated_groups["Junctions"][n_groups] = junctions_subset_new_group
            identified_parallels.extend([parallels_j, parallels_i] + junctions_subset_new_group)
            n_groups += 1
    else : # if p_i has no crossroad in any other parallels
        associated_groups["Main_roads"][n_groups] = [parallels_i]
        if parallels_i in parallels_with_junctions.keys():
            junctions_subset_new_group = [parallels_j for parallels_j in parallels_with_junctions[parallels_i] if parallels_j not in identified_parallels]
            associated_groups["Junctions"][n_groups] = junctions_subset_new_group
            identified_parallels.extend(junctions_subset_new_group)
        else :
            associated_groups["Junctions"][n_groups] = []
        n_groups += 1
    identified_parallels.append(parallels_i)
groups_ids = list(np.arange(0,n_groups))
time_total = time.time() - start_3
time_t = time.gmtime(time_total)
line = "It took" +  str(time_t.tm_mday )+ "days"+ str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

line="3.b) Create the groups and get the road ids within their subsets of Main_roads and Junctions and save in Groups and Roads_to_Groups"
print(line)
f.write(line + '\n')
# NB : For now, subsets of main roads and junctions are identified with their parallels id
Groups = {}
for kth_group in range(n_groups):
    Groups[kth_group] = {'MAIN_ROAD_IDS': [], 'JUNCTION_IDS': [], 'NODES_IDS': [], 'JUNCTION_NODES': [], 'TUNNEL':False}
    group_main_roads = []
    group_junctions = []
    group_nodes_ids = []
    junction_nodes_ids = []
    group_max_length = 0
    # get group information from previously identified parallels
    for parallels_p in associated_groups["Main_roads"][kth_group] :
        group_main_roads.extend(roads_of_parallels_id[parallels_p])
        group_nodes_ids.extend(parallels_nodes[parallels_p])
        if parallels_max_length[parallels_p] > group_max_length :
            group_max_length = parallels_max_length[parallels_p]
    for parallels_c in associated_groups["Junctions"][kth_group] :
        group_junctions.extend(roads_of_parallels_id[parallels_c])
        junction_nodes_ids.extend(parallels_nodes[parallels_c])
        if parallels_max_length[parallels_c] > group_max_length :
            group_max_length = parallels_max_length[parallels_c]
    # save information in Groups and Roads_to_groups
    group_roads = group_junctions + group_main_roads
    Groups[kth_group]["MAIN_ROAD_IDS"] = group_main_roads
    Groups[kth_group]["JUNCTION_IDS"] = group_junctions
    Groups[kth_group]["NODES_IDS"] = list(set(group_nodes_ids))
    Groups[kth_group]["JUNCTION_NODES"] = list(set(junction_nodes_ids))
    Groups[kth_group]["TUNNEL"] = tunnel[parallels_p]
    for road_i in group_roads :
        i = roads_ids.index(road_i)
        Roads_to_groups.loc[i, 'GROUP_ID'] = kth_group
        Roads_to_groups.loc[i, 'COLOR'] = kth_group % number_of_colors
        Roads_to_groups.loc[i, 'MERGED'] = True
        if road_i in group_main_roads :
            Roads_to_groups.loc[i, 'MAIN_ROAD'] = True
        else :
            Roads_to_groups.loc[i, 'MAIN_ROAD'] = False
time_total = time.time() - start_3
time_t = time.gmtime(time_total)
line = "It took" +  str(time_t.tm_mday )+ "days"+ str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

line="3.c) Check if any non merged road is a junction of a main road"
print(line)
f.write(line + '\n')
non_merged_roads = non_merged_roads + roundabouts_ids # add roundabouts
for cell_id in list_of_cells :
    if len(parallels_per_cell[cell_id]) == 0:
        continue
    potential_junctions = [item for item in roads_per_cell[cell_id] if (item in non_merged_roads)]
    for road_i in potential_junctions :
        non_merged_roads.remove(road_i)
        i = roads_ids.index(road_i)
        node_A = str(roads.loc[i, 'NODE_A'])
        node_B = str(roads.loc[i, 'NODE_B'])
        # filter groups that are close enough
        neighborhood_links_i = neighborhood_links[road_i]
        potential_groups_i = Roads_to_groups[(Roads_to_groups.ROAD_ID.isin(neighborhood_links_i)) & (Roads_to_groups.MAIN_ROAD == True)].reset_index(drop=True)
        group_number_max = potential_groups_i.shape[0]
        merged = False
        if group_number_max > 0 :
            groups_cnt = 0
            # for each group to merge, check if road_i intersects at least two of the group's roads
            while not (merged or groups_cnt == group_number_max):
                # check if road_i is a crossroad of this road
                group_number = potential_groups_i.loc[groups_cnt, 'GROUP_ID']
                nodes = [node_A, node_B] # list of nodes that have not been identified as connected to the group yet
                for node in Groups[group_number]["NODES_IDS"] :
                    if (str(node) == node_A) and (node_A in nodes) :
                        nodes.remove(node_A)
                    elif (str(node) == node_B) and (node_B in nodes) :
                        nodes.remove(node_B)
                if len(nodes) == 0 : # i.e. if both nodes are connected to the group's parallel roads
                    # update the group's roads
                    group_i_junctions = Groups[group_number]["JUNCTION_IDS"]
                    group_i_junctions.append(road_i)
                    Groups[group_number]["JUNCTION_IDS"] = group_i_junctions
                    # update Roads_to_groups
                    Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[i, 'MAIN_ROAD'] = False
                    Roads_to_groups.loc[i, 'MERGED'] = True
                    merged = True
                else :
                    groups_cnt = groups_cnt + 1
        if not merged :
            non_merged_roads.append(road_i)
time_total = time.time() - start_3
time_t = time.gmtime(time_total)
line = "It took" +  str(time_t.tm_mday )+ "days"+ str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

# 3.d) Ajout des groupes de cardinal 1
line="3.d) For identified streets strictly equal to a road, save their information"
print(line)
f.write(line + '\n')
n_street_roads = len(non_merged_roads)
print("number of groups ids,", n_groups, ", number of street-roads", n_street_roads)
for k in range(n_street_roads):
    link_k = non_merged_roads[k]
    line = str(link_k) + ', type:' + str(type(link_k))
    f.write(line + '\n')
    road_index = Roads_to_groups[Roads_to_groups.ROAD_ID == str(link_k)].index[0]
    kth_group = n_groups + 1 + k
    # update Roads_to_groups
    Roads_to_groups.loc[road_index, 'GROUP_ID'] = kth_group
    Roads_to_groups.loc[road_index, 'MAIN_ROAD'] = True
    Roads_to_groups.loc[road_index, 'MERGED'] = False
    Roads_to_groups.loc[road_index, 'COLOR'] = kth_group % number_of_colors
    node_A = roads.loc[road_index,'NODE_A']
    node_B = roads.loc[road_index,'NODE_B']
    # save information in Groups
    Groups[kth_group] = {'MAIN_ROAD_IDS': [], 'CROSSROAD_IDS': [], 'NODES_IDS': [], 'CROSSROAD_NODES': [], 'TUNNEL':False}
    Groups[kth_group]['NODES_IDS'] = [node_A, node_B]
    Groups[kth_group]['JUNCTION_IDS'] = []
    Groups[kth_group]['JUNCTION_NODES'] = []
    Groups[kth_group]['MAIN_ROAD_IDS'] = [link_k]
    Groups[kth_group]['TUNNEL'] = Roads_to_groups.loc[road_index, 'TUNNEL']
time_total = time.time() - start_3
time_t = time.gmtime(time_total)
line = "It took" + str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

# 4. Get groups information for nodes (4h)
line="4. Get groups information for nodes"
print(line)
f.write(line + '\n')
Nodes = {'NODE_ID':[], 'X':[], 'Y':[], 'ROAD_LINKS':[],'GROUPS':[], 'N_GROUPS':[]}
for i, node_i_links in enumerate(nodes_OSM['LINKS']):
    node_i_groups = []
    for link in node_i_links :
        if str(link) in roads_ids :
            link_index = Roads_to_groups[Roads_to_groups.ROAD_ID == str(link)].index[0]
            new_group = Roads_to_groups.loc[link_index, 'GROUP_ID']
            node_i_groups.append(new_group)
    node_i_groups = list(set(node_i_groups))
    if len(node_i_groups) > 0:
        Nodes['NODE_ID'].append(nodes_OSM.loc[i,'ID'])
        Nodes['X'].append(nodes_OSM.loc[i,'X'])
        Nodes['Y'].append(nodes_OSM.loc[i,'Y'])
        Nodes['ROAD_LINKS'].append(node_i_links)
        Nodes['GROUPS'].append(node_i_groups)
        Nodes['N_GROUPS'].append(len(node_i_groups))
    if len(node_i_groups) == 0:
        print(i, node_i_links)
        print(nodes_OSM['ID'][i])

time_total = time.time() - start_parallel
time_t = time.gmtime(time_total)
line = "Total time : It took" + str(time_t.tm_hour)+ 'hour'+ str(time_t.tm_min)+ 'min'+ str(time_t.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

# Save direct and inverse link transform
n_streets = n_street_roads + n_groups
print('number of non merged roads',n_street_roads)
print("Number of streets: ",n_streets)
print("Number of groups: ",n_groups)
pickle.dump(Groups, open(fout_groups, 'wb'))
Roads_to_groups.to_csv(fout_roads, index=False)
pickle.dump(Roads_to_groups, open(fout_roads_dat, 'wb'))
pickle.dump(Nodes, open(fout_nodes, 'wb'))
f.close()


######################
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Polygon,Point,LineString
import numpy as np
from include.create_poly_around_road import create_1_poly_around_road
from include.create_circles import create_circles_around_road_nodes
from include.lines_geometry import *
from include.get_indexes import get_indexes
import pickle,sys
import time

# Load grid data
path = "../../"
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
grid=pickle.load(open(gridfile,'rb'))

# Get cell information
poly_num = 4672
x0 = grid['we_id'][poly_num]
y0 = grid['sn_id'][poly_num]
prefix = 'we_id_from_' + str(x0) + '_sn_id_from_' + str(y0)

# Input and output files
dir_roads = path + "temp/load_OSM/subdomains_roads_lamb93/" # il faut prendre le outdir de road_graph
detailed_roads_f = path + "temp/load_OSM/roads_links.csv"
nodes_f = path + "temp/load_OSM/roads_nodes.csv"
# buildings_f = path + 'temp/data/load_BDTOPO/buildings.dat'
dir_buildings = path + 'temp/load_BDTOPO/subdomains_buildings_with_buffer_100m_lamb93/'
waterways_f = path + 'temp/load_OSM/waterways/filtered_waterways.shp'
fname = 'sub_' + prefix + '.dat'
outdir_groups = path + 'temp/street_graph/groups/'
fout_nodes = outdir_groups + 'group_nodes.dat'
fout_roads = outdir_groups + 'roads_to_groups.csv'
fout_groups = outdir_groups + 'groups.dat'

# Load data
roads = pickle.load(open(dir_roads + fname, 'rb'))
roads['ROAD_ID'] = [str(road) for road in roads['ROAD_ID']]
roads_df = pd.DataFrame(roads)
detailed_roads = pd.read_csv(detailed_roads_f)
detailed_roads = detailed_roads.astype({'ID':'str'})
roads_ids = list(roads['ROAD_ID'])
#roads_df = pd.merge(left = roads_df, right=detailed_roads, how='left', on=['ROAD_ID'])
nodes_types = {'ID':str,'X':float, 'Y':float, 'LINKS':object}
nodes_OSM = pd.read_csv(nodes_f, dtype=nodes_types)
print(nodes_OSM['LINKS'])
nodes_OSM['LINKS'] = nodes_OSM['LINKS'].apply(eval)
buildings_data = pickle.load(open(dir_buildings + fname, 'rb'))
buildings = gpd.GeoDataFrame(data=buildings_data, geometry=buildings_data['POLYGON'])
waterways = gpd.read_file(waterways_f)
print('Data loaded... ', time.time())

# Initialize data
Roads_to_groups = pd.DataFrame(columns=['ROAD_ID', 'GROUP_ID', 'PARALLEL', 'MERGED', 'COLOR', 'ROUNDABOUT', 'TUNNEL', 'GEOMETRY'])
# Groups = {'GROUP_ID':[], 'ROAD_IDS':[], 'PARALLEL_ROADS_IDS':[], 'NODES_IDS':[], 'MAX_LENGTH':[]}
Groups = {}
roads_list = list(detailed_roads['ID'])
for i, road_i in enumerate(roads_ids) :
    road_index = roads_list.index(road_i)
    Roads_to_groups.loc[i,'ROAD_ID'] = road_i
    Roads_to_groups.loc[i,'ROUNDABOUT'] = detailed_roads.loc[road_index,'ROUNDABOUT']
    Roads_to_groups.loc[i,'TUNNEL'] = detailed_roads.loc[road_index,'TUNNEL']
    Roads_to_groups.loc[i,'GEOMETRY'] = detailed_roads.loc[road_index,'GEOMETRY']

# Parameters: thresholds for merging
min_radius_same_points = 5 # meters
poly_width = 100  # meters
building_area_threshold = 0.05
road_length_percentage = 0.8
number_of_colors = 25
building_height_threshold = 5 # meters
# small_group_length = 40 # meters

# Filter round-abouts
Roundabouts = Roads_to_groups[Roads_to_groups.ROUNDABOUT == True]
# Roads_to_groups = Roads_to_groups[Roads_to_groups.ROUNDABOUT == False]
roundabouts_ids = list(Roundabouts['ROAD_ID'])

# 1. Identify parallel roads
print('1. Identify parallel roads')
start_parallel = time.time()
merging_roads_groups = {} # output: dictionary of groups of roads
parallel_roads = {}
group_nodes = {}
longest_road_within_group = {}
groups_max_length = {}
tunnel = {}
roads_ids_without_roundabouts = [item for item in roads_ids if item not in roundabouts_ids]
non_merged_roads = roads_ids_without_roundabouts.copy() # this list will be modified
for road_i in roads_ids_without_roundabouts :
    i = roads_ids.index(road_i)
    # purpose: to test which roads have to be added to the group of roads to merge with road_i, if there is
    if road_i in non_merged_roads :
        # Get road geometry
        node_A = roads['NODE_A'][i]
        node_B = roads['NODE_B'][i]
        road_xa = roads['XA'][i]
        road_ya = roads['YA'][i]
        road_xb = roads['XB'][i]
        road_yb = roads['YB'][i]
        points_coordinates_1 = [road_xa, road_ya, road_xb, road_yb]
        length_i = roads["LENGTH"][i]
        tunnel_i = Roads_to_groups.loc[i, 'TUNNEL']
        # Define circles around road_i nodes and a polygon.
        radius = min(min_radius_same_points, length_i)
        circle_A, circle_B = create_circles_around_road_nodes(road_xa, road_xb, road_ya, road_yb, radius)
        poly_i = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, poly_width)
        # Initialize merging
        non_merged_roads.remove(road_i)
        merged = False
        n_groups = len(longest_road_within_group)
    else :
        merged = True
    roads_to_test = list(longest_road_within_group.values()) + non_merged_roads  # add all roads that are not merged for the moment
    cnt = 0 # counter of the place within roads_to_test
    while not (merged or cnt==len(roads_to_test)):
        road_j = roads_to_test[cnt]
        j = roads_ids.index(road_j)
        tunnel_j = Roads_to_groups.loc[j, 'TUNNEL']
        if tunnel_i == tunnel_j :
            # data for tests
            road_xc = roads['XA'][j]
            road_yc = roads['YA'][j]
            road_xd = roads['XB'][j]
            road_yd = roads['YB'][j]
            node_C_point = Point(road_xc, road_yc)
            node_D_point = Point(road_xd, road_yd)
            node_C = roads['NODE_A'][j]
            node_D = roads['NODE_B'][j]
            edge_j = roads["GEOMETRY"][j]
            length_j = roads["LENGTH"][j]
            points_coordinates_2 = [road_xc, road_yc, road_xd, road_yd]
            if poly_i.intersects(edge_j):
                intersect_i_j = poly_i.intersection(edge_j).length
            else:
                intersect_i_j = 0
            ## Test 1: create_circles
            if node_C_point.intersects(circle_A) and node_D_point.intersects(circle_B) and parallel_lines(points_coordinates_1, points_coordinates_2):
                # adding road_i in merging_roads_groups to the corresponding list
                if (n_groups > 0) and (cnt + 1 <= n_groups) : # road_j is already in a group, add road_i to it
                    # cnt is the index of the longest_road_within_group
                    group_number = cnt
                    merging_roads_groups[group_number].append(road_i)
                    parallel_roads[group_number].append(road_i)
                    group_nodes[group_number].extend([node_A, node_B])
                    Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[i, 'PARALLEL'] = True
                    Roads_to_groups.loc[i, 'MERGED'] = True
                    Roads_to_groups.loc[i, 'STEP'] = 1
                    if length_i > length_j :
                        longest_road_within_group[group_number] = road_i
                        groups_max_length[group_number] = length_i
                else : # else, create a new group with road_i and road_j
                    group_number = n_groups
                    merging_roads_groups[group_number] = [road_i, road_j]
                    parallel_roads[group_number] = [road_i, road_j]
                    group_nodes[group_number] = [node_A, node_B, node_C, node_D]
                    tunnel[group_number] = tunnel_i
                    Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[i, 'PARALLEL'] = True
                    Roads_to_groups.loc[i, 'MERGED'] = True
                    Roads_to_groups.loc[i, 'STEP'] = 1.1
                    Roads_to_groups.loc[j, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[j, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[j, 'PARALLEL'] = True
                    Roads_to_groups.loc[j, 'MERGED'] = True
                    Roads_to_groups.loc[j, 'STEP'] = 1
                    non_merged_roads.remove(road_j)
                    if length_j > length_i :
                        longest_road_within_group[group_number] = road_j
                        groups_max_length[group_number] = length_j
                    else :
                        longest_road_within_group[group_number] = road_i
                        groups_max_length[group_number] = length_i
                merged = True
            elif node_C_point.intersects(circle_B) and node_D_point.intersects(circle_A) and parallel_lines(points_coordinates_1, points_coordinates_2):
                # adding road_i in merging_roads_groups to the corresponding list
                if (n_groups > 0) and (cnt + 1 <= n_groups):  # road_j is already in a group, add road_i to it
                    group_number = cnt
                    merging_roads_groups[group_number].append(road_i)
                    parallel_roads[group_number].append(road_i)
                    group_nodes[group_number].extend([node_A, node_B])
                    Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[i, 'PARALLEL'] = True
                    Roads_to_groups.loc[i, 'MERGED'] = True
                    Roads_to_groups.loc[i, 'STEP'] = 1
                    if length_i > length_j :
                        longest_road_within_group[group_number] = road_i
                        groups_max_length[group_number] = length_i
                else:  # else, create a new group with road_i and road_j
                    group_number = n_groups
                    merging_roads_groups[group_number] = [road_i, road_j]
                    parallel_roads[group_number] = [road_i, road_j]
                    group_nodes[group_number] = [node_A, node_B, node_C, node_D]
                    tunnel[group_number] = tunnel_i
                    Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[i, 'PARALLEL'] = True
                    Roads_to_groups.loc[i, 'MERGED'] = True
                    Roads_to_groups.loc[i, 'STEP'] = 1
                    Roads_to_groups.loc[j, 'GROUP_ID'] = group_number
                    Roads_to_groups.loc[j, 'COLOR'] = group_number % number_of_colors
                    Roads_to_groups.loc[j, 'PARALLEL'] = True
                    Roads_to_groups.loc[j, 'MERGED'] = True
                    Roads_to_groups.loc[j, 'STEP'] = 1
                    non_merged_roads.remove(road_j)
                    if length_j > length_i:
                        longest_road_within_group[group_number] = road_j
                        groups_max_length[group_number] = length_j
                    else:
                        longest_road_within_group[group_number] = road_i
                        groups_max_length[group_number] = length_i
                merged = True
            ## Test 2: poly_around_road (under the condition: no common node)
            elif ((intersect_i_j > road_length_percentage * length_i) or (intersect_i_j > road_length_percentage * length_j)) and (node_C not in [node_A, node_B]) and (node_D not in [node_A, node_B]) and parallel_lines(points_coordinates_1, points_coordinates_2):
                # Check that there is no building between the two roads (condition on area)
                if shapely.distance(Point(road_xa,road_ya), Point(road_xc,road_yc)) > shapely.distance(Point(road_xa,road_ya), Point(road_xd,road_yd)) :
                    poly_roads = Polygon(((road_xa, road_ya),(road_xb, road_yb),(road_xc, road_yc),(road_xd, road_yd)))
                else :
                    poly_roads = Polygon(((road_xa, road_ya), (road_xb, road_yb), (road_xd, road_yd), (road_xc, road_yc)))
                subdata = buildings[buildings.geometry.intersects(poly_roads)]
                subdata = subdata[subdata.HEIGHT > building_height_threshold]
                total_buildings_area = sum(subdata.geometry.area)
                intersecting_waterways = waterways[waterways.geometry.intersects(poly_roads)]
                # previous option: len(list(subdata.geometry)) == 0
                if (total_buildings_area < building_area_threshold * poly_roads.area) and (
                        sum(intersecting_waterways.geometry.length) == 0): # if there is no building large enough, then merge the two roads
                    # adding road_i in merging_roads_groups to the corresponding list
                    if (n_groups > 0) and (cnt + 1 <= n_groups):  # road_j is already in a group, add road_i to it
                        group_number = cnt
                        merging_roads_groups[group_number].append(road_i)
                        parallel_roads[group_number].append(road_i)
                        group_nodes[group_number].extend([node_A, node_B])
                        Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                        Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                        Roads_to_groups.loc[i, 'PARALLEL'] = True
                        Roads_to_groups.loc[i, 'MERGED'] = True
                        Roads_to_groups.loc[i, 'STEP'] = 1
                        if length_i > length_j :
                            longest_road_within_group[group_number] = road_i
                            groups_max_length[group_number] = length_i
                    else:  # else, create a new group with road_i and road_j
                        group_number = n_groups
                        merging_roads_groups[group_number] = [road_i, road_j]
                        parallel_roads[group_number] = [road_i, road_j]
                        group_nodes[group_number] = [node_A, node_B, node_C, node_D]
                        tunnel[group_number] = tunnel_i
                        Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                        Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                        Roads_to_groups.loc[i, 'PARALLEL'] = True
                        Roads_to_groups.loc[i, 'MERGED'] = True
                        Roads_to_groups.loc[i, 'STEP'] = 1
                        Roads_to_groups.loc[j, 'GROUP_ID'] = group_number
                        Roads_to_groups.loc[j, 'COLOR'] = group_number % number_of_colors
                        Roads_to_groups.loc[j, 'PARALLEL'] = True
                        Roads_to_groups.loc[j, 'MERGED'] = True
                        Roads_to_groups.loc[j, 'STEP'] = 1
                        non_merged_roads.remove(road_j)
                        if length_j > length_i:
                            longest_road_within_group[group_number] = road_j
                            groups_max_length[group_number] = length_j
                        else:
                            longest_road_within_group[group_number] = road_i
                            groups_max_length[group_number] = length_i
                    merged = True
                else :
                    cnt = cnt + 1
            else :
                cnt = cnt + 1
        else :
            cnt = cnt + 1
    if not merged :
        non_merged_roads.append(road_i)
        Roads_to_groups.loc[i, 'MERGED'] = False
n_groups = len(merging_roads_groups)
groups_ids = list(np.arange(0,n_groups))

time_parallell_roads_sec = time.time()-start_parallel
time_PR = time.gmtime(time_parallell_roads_sec)
print("It took", time_PR.tm_hour, 'hour', time_PR.tm_min, 'min', time_PR.tm_sec ,'sec.')
# In the end : non_merged_roads and merging_roads_groups with groups of roads

# 2. Check if 2 groups are to be merged, based on longest road of each group
groups_rep_road = list(longest_road_within_group.values())
roads_to_test = list(longest_road_within_group.values())
for group_i, road_i in enumerate(groups_rep_road) :
    i = roads_ids.index(road_i)
    # purpose: to test which roads have to be added to the group of roads to merge with road_i, if there is
    if road_i in roads_to_test :
        # Get road geometry
        node_A = roads['NODE_A'][i]
        node_B = roads['NODE_B'][i]
        road_xa = roads['XA'][i]
        road_ya = roads['YA'][i]
        road_xb = roads['XB'][i]
        road_yb = roads['YB'][i]
        points_coordinates_1 = [road_xa, road_ya, road_xb, road_yb]
        length_i = roads["LENGTH"][i]
        tunnel_i = tunnel[group_i]# Roads_to_groups.loc[i, 'TUNNEL']
        # Define circles around road_i nodes and a polygon.
        radius = min(min_radius_same_points, length_i)
        circle_A, circle_B = create_circles_around_road_nodes(road_xa, road_xb, road_ya, road_yb, radius)
        poly_i = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, poly_width)
        # Initialize merging
        merged = False
    else :
        merged = True
    group_j = 0 # counter of the place within roads_to_test
    while not (merged or group_j == n_groups): # len(roads_to_test)
        if (group_j == group_i) or (group_j not in groups_ids) : # or (group_j not in roads_to_test)
            group_j = group_j + 1
            continue
        tunnel_j = tunnel[group_j]
        if tunnel_i != tunnel_j :
            group_j = group_j + 1
            continue
        road_j = longest_road_within_group[group_j]  # longest
        j = roads_ids.index(road_j)
        # data for tests
        road_xc = roads['XA'][j]
        road_yc = roads['YA'][j]
        road_xd = roads['XB'][j]
        road_yd = roads['YB'][j]
        node_C_point = Point(road_xc, road_yc)
        node_D_point = Point(road_xd, road_yd)
        node_C = roads['NODE_A'][j]
        node_D = roads['NODE_B'][j]
        edge_j = roads["GEOMETRY"][j]
        length_j = roads["LENGTH"][j]
        points_coordinates_2 = [road_xc, road_yc, road_xd, road_yd]
        if poly_i.intersects(edge_j):
            intersect_i_j = poly_i.intersection(edge_j).length
        else:
            intersect_i_j = 0
        # test road_i vs road_j : circle, then rectangles....
        ## Test 1: create_circles
        if node_C_point.intersects(circle_A) and node_D_point.intersects(circle_B) and parallel_lines(
                points_coordinates_1, points_coordinates_2): # then we only keep group_j and merge it with group_i
            if length_i > length_j:
                longest_road_within_group[group_j] = road_i
                groups_max_length[group_j] = length_i
                roads_to_test.remove(road_j)
            else :
                roads_to_test.remove(road_i)
            merged = True
            parallel_roads[group_j] = list(set(parallel_roads[group_i] + parallel_roads[group_j]))
            group_nodes[group_j] = list(set(group_nodes[group_i] + group_nodes[group_j]))
            merging_roads_groups[group_j] = list(set(merging_roads_groups[group_i] + merging_roads_groups[group_j]))
            for road_k in merging_roads_groups[group_i] :
                k = roads_ids.index(road_k)
                Roads_to_groups.loc[k, 'GROUP_ID'] = group_j
                Roads_to_groups.loc[k, 'COLOR'] = group_j % number_of_colors
            parallel_roads.pop(group_i)
            group_nodes.pop(group_i)
            merging_roads_groups.pop(group_i)
            groups_ids.remove(group_i)
            groups_max_length.pop(group_i)
            longest_road_within_group.pop(group_i)
            tunnel.pop(group_i)
        elif node_C_point.intersects(circle_B) and node_D_point.intersects(circle_A) and parallel_lines(
            points_coordinates_1, points_coordinates_2): # then we only keep group_j and merge it with group_i
            if length_i > length_j:
                longest_road_within_group[group_j] = road_i
                groups_max_length[group_j] = length_i
                roads_to_test.remove(road_j)
            else:
                roads_to_test.remove(road_i)
            merged = True
            parallel_roads[group_j] = list(set(parallel_roads[group_i] + parallel_roads[group_j]))
            group_nodes[group_j] = list(set(group_nodes[group_i] + group_nodes[group_j]))
            merging_roads_groups[group_j] = list(set(merging_roads_groups[group_i] + merging_roads_groups[group_j]))
            for road_k in merging_roads_groups[group_i]:
                k = roads_ids.index(road_k)
                Roads_to_groups.loc[k, 'GROUP_ID'] = group_j
                Roads_to_groups.loc[k, 'COLOR'] = group_j % number_of_colors
            parallel_roads.pop(group_i)
            group_nodes.pop(group_i)
            merging_roads_groups.pop(group_i)
            groups_ids.remove(group_i)
            groups_max_length.pop(group_i)
            longest_road_within_group.pop(group_i)
            tunnel.pop(group_i)
        ## Test 2: poly_around_road (under the condition: no common node)
        elif ((intersect_i_j > road_length_percentage * length_i) or (intersect_i_j > road_length_percentage * length_j)) and (node_C not in [node_A, node_B]) and (
            node_D not in [node_A, node_B]) and parallel_lines(points_coordinates_1, points_coordinates_2): # then we only keep group_j and merge it with group_i
            # Check that there is no building between the two roads (condition on area)
            if shapely.distance(Point(road_xa, road_ya), Point(road_xc, road_yc)) > shapely.distance(
                    Point(road_xa, road_ya), Point(road_xd, road_yd)):
                poly_roads = Polygon(((road_xa, road_ya), (road_xb, road_yb), (road_xc, road_yc), (road_xd, road_yd)))
            else:
                poly_roads = Polygon(((road_xa, road_ya), (road_xb, road_yb), (road_xd, road_yd), (road_xc, road_yc)))
            subdata = buildings[buildings.geometry.intersects(poly_roads)]
            subdata = subdata[subdata.HEIGHT > building_height_threshold]
            total_buildings_area = sum(subdata.geometry.area)
            intersecting_waterways = waterways[waterways.geometry.intersects(poly_roads)]
            # previous option: len(list(subdata.geometry)) == 0
            if (total_buildings_area < building_area_threshold * poly_roads.area)  and (sum(intersecting_waterways.geometry.length) == 0):  # if there is no building large enough, then merge the two roads
                # adding road_i in merging_roads_groups to the corresponding list
                if length_i > length_j:
                    longest_road_within_group[group_j] = road_i
                    groups_max_length[group_j] = length_i
                    roads_to_test.remove(road_j)
                else:
                    roads_to_test.remove(road_i)
                merged = True
                parallel_roads[group_j] = list(set(parallel_roads[group_i] + parallel_roads[group_j]))
                group_nodes[group_j] = list(set(group_nodes[group_i] + group_nodes[group_j]))
                merging_roads_groups[group_j] = list(set(merging_roads_groups[group_i] + merging_roads_groups[group_j]))
                for road_k in merging_roads_groups[group_i]:
                    k = roads_ids.index(road_k)
                    Roads_to_groups.loc[k, 'GROUP_ID'] = group_j
                    Roads_to_groups.loc[k, 'COLOR'] = group_j % number_of_colors
                parallel_roads.pop(group_i)
                group_nodes.pop(group_i)
                merging_roads_groups.pop(group_i)
                groups_ids.remove(group_i)
                groups_max_length.pop(group_i)
                longest_road_within_group.pop(group_i)
                tunnel.pop(group_i)
            else :
                group_j = group_j + 1
        else:
            group_j = group_j + 1


########################## Perpendicular groups #################

print("3. Perpendicular groups")
groups_crossroads = pd.DataFrame(columns = ['INITIAL_GROUP_ID', 'GROUP_ID', 'LONGEST'])
print("groups_ids", groups_ids)
test_merging_roads_groups = merging_roads_groups.copy()
## First, for each road of a group, check if it is a crossroad
for group_i, roads_within_group in test_merging_roads_groups.items() : # small_merging_roads_groups = list of groups, each of them is a list of the group's road_ids
    groups_to_check = groups_ids.copy()
    for road_i in roads_within_group :
        i = roads_ids.index(road_i)
        node_A = str(roads['NODE_A'][i])
        node_B = str(roads['NODE_B'][i])
        merged = False
        group_j = 0
        # for each group to merge, check if road_i intersects at least two of the group's roads
        while not (merged or group_j == n_groups):
            if (group_j in groups_to_check) and (group_i != group_j) and (tunnel[group_i] == tunnel[group_j]) :
                nodes = [node_A, node_B]  # list of nodes that have not been identified as connected to the group yet
                for node in group_nodes[group_j]:
                    if (str(node) == node_A) and (node_A in nodes):
                        nodes.remove(node_A)
                    elif (str(node) == node_B) and (node_B in nodes):
                        nodes.remove(node_B)
                if len(nodes) == 0:  # i.e. if both nodes are connected to the group's parallel roads
                    merged = True
                    # if group_j has the longest road, keep it
                    if groups_max_length[group_i] < groups_max_length[group_j] :
                        # groups_i_j = [[road_i, group_i, group_j]]
                        groups_i_j = [[group_i, group_j, group_i]]
                        row = pd.DataFrame(data = groups_i_j, columns=['INITIAL_GROUP_ID', 'GROUP_ID', 'LONGEST'])
                        groups_crossroads = pd.concat([groups_crossroads, row], ignore_index=True)
                        groups_to_check.remove(group_j)
                    else : # else keep group_i
                        # road_j = groups_i_j = [[road_j, group_j, group_i]]
                        road_j = longest_road_within_group[group_j]
                        groups_i_j = [[group_i, group_j, group_j]]
                        row = pd.DataFrame(data=groups_i_j, columns=['INITIAL_GROUP_ID', 'GROUP_ID', 'LONGEST'])
                        groups_crossroads = pd.concat([groups_crossroads, row], ignore_index=True)
                        groups_to_check.remove(group_j)
                else:
                    group_j += 1
            else :
                group_j += 1
groups_crossroads_list = groups_crossroads['INITIAL_GROUP_ID']
groups_crossroads_identification = {}
for g in groups_crossroads_list :
    g_crossroad_of_groups = groups_crossroads[groups_crossroads.INITIAL_GROUP_ID == g]
    groups_crossroads_identification[g] = list(g_crossroad_of_groups['GROUP_ID'])

# identify groups such as g1 has a crossroad in g2 and vice-versa (Property 2.2)
perpendicularity = {}
for g1, g2_list in groups_crossroads_identification.items() :
    for g2 in g2_list :
        if g2 in groups_crossroads_identification.keys() and g1 in groups_crossroads_identification[g2] : # i.e. g1 perp. g2
            # check which one is the longest
            if groups_max_length[g1] < groups_max_length[g2] :
                if g1 in perpendicularity.keys() :
                    if g2 not in perpendicularity[g1] :
                        perpendicularity[g1].append(g2)
                        groups_crossroads_identification[g2].remove(g1)
                else :
                    perpendicularity[g1] = [g2]
                    groups_crossroads_identification[g2].remove(g1)
            else :
                if g2 in perpendicularity.keys() :
                    if g1 not in perpendicularity[g2] :
                        perpendicularity[g2].append(g1)
                        groups_crossroads_identification[g1].remove(g2)
                else:
                    perpendicularity[g2] = [g1]
                    groups_crossroads_identification[g1].remove(g2)

# Deduce perpendicularity for all groups (Theorem 2.1)
for g1 in groups_crossroads_identification.keys():
    if g1 in perpendicularity.keys() :
        perpendicularity[g1] = perpendicularity[g1][0]
    else :
        if len(groups_crossroads_identification[g1]) > 0 :
            l = 0
            for g2 in groups_crossroads_identification[g1] :
                if groups_max_length[g2] > l :
                    l = groups_max_length[g2]
                    perpendicularity[g1] = g2

# Identify parallel groups
parallelism = {}
perpendicularity_initial = perpendicularity.copy()
for g1, g2 in perpendicularity_initial.items() : # if g1 perp. g2
    if g2 in perpendicularity_initial.keys() : # and g2 perp. g3
        g3 = perpendicularity_initial[g2]
        parallelism[g1] = g3 # then g1 perp. g3
        perpendicularity.pop(g1)

def check_perpendicular_groups(perpendicular:dict, parallel:dict) :
    # Identify perpendicular groups to parallel groups
    for g1, g2 in perpendicular.items():  # if g1 perp. g2
        if g2 in parallel.keys():  # and g2 parallel g3
            g3 = parallel[g2]
            perpendicular[g1] = g3  # then g1 perp. g3
    # Identify parallel groups to perpendicular groups
    parallelism_initial = parallel.copy()
    for g1, g2 in parallelism_initial.items():  # if g1 perp. g2
        if g2 in perpendicular.keys():  # and g2 parallel g3
            g3 = perpendicular[g2]
            perpendicular[g1] = g3  # then g1 perp. g3
            parallel.pop(g1)
    return(perpendicular, parallel)

modification_1 = any(check in parallelism.keys() for check in perpendicularity.values())
modification_2 = any(check in perpendicularity.keys() for check in parallelism.values())
while (modification_1 or modification_2) :
    perpendicularity, parallelism = check_perpendicular_groups(perpendicularity, parallelism)
    modification_1 = any(check in parallelism.keys() for check in perpendicularity.values())
    modification_2 = any(check in perpendicularity.keys() for check in parallelism.values())

print("perpendicularity",perpendicularity)
print("parallel",parallelism)

## Then for remaining roads, check if the road(s) they were associated with is a crossroad. If it is, then they both are.
for initial_group_number, group_i in perpendicularity.items() :
    merging_roads_groups[group_i] = list(set(merging_roads_groups[group_i] + merging_roads_groups[initial_group_number]))
    group_nodes[group_i] = list(set(group_nodes[group_i] + group_nodes[initial_group_number]))
    # parallel_roads[group_i]  already complete
    for road_k in merging_roads_groups[initial_group_number] :
        i = roads_ids.index(road_k)
        Roads_to_groups.loc[i, 'GROUP_ID'] = group_i
        Roads_to_groups.loc[i, 'COLOR'] = group_i % number_of_colors
        Roads_to_groups.loc[i, 'PARALLEL'] = False
        Roads_to_groups.loc[i, 'MERGED'] = True
        Roads_to_groups.loc[i, 'STEP'] = 2
    merging_roads_groups.pop(initial_group_number)
    groups_ids.remove(initial_group_number)
    parallel_roads.pop(initial_group_number)
    groups_max_length.pop(initial_group_number)
    longest_road_within_group.pop(initial_group_number)

#  check if the road(s) they were associated with is in a parallel groups
for initial_group_number, group_i in parallelism.items() :
    merging_roads_groups[group_i] = list(set(merging_roads_groups[group_i] + merging_roads_groups[initial_group_number]))
    group_nodes[group_i] = list(set(group_nodes[group_i] + group_nodes[initial_group_number]))
    parallel_roads[group_i] = list(set(parallel_roads[group_i] + parallel_roads[initial_group_number]))
    for road_k in merging_roads_groups[initial_group_number] :
        i = roads_ids.index(road_k)
        Roads_to_groups.loc[i, 'GROUP_ID'] = group_i
        Roads_to_groups.loc[i, 'COLOR'] = group_i % number_of_colors
        Roads_to_groups.loc[i, 'PARALLEL'] = True
        Roads_to_groups.loc[i, 'MERGED'] = True
        Roads_to_groups.loc[i, 'STEP'] = 2
    merging_roads_groups.pop(initial_group_number)
    groups_ids.remove(initial_group_number)
    parallel_roads.pop(initial_group_number)
    groups_max_length.pop(initial_group_number)
    longest_road_within_group.pop(initial_group_number)
    tunnel.pop(initial_group_number)
# update groups ids and the number of groups
groups_ids = list(groups_max_length.keys())
n_groups = len(groups_ids)
print("groups_ids", groups_ids)

########################## Crossroads identification #################

print("4. Identify crossroads")
start = time.time()
non_merged_roads = non_merged_roads + roundabouts_ids # add roundabouts
group_number_max = max(groups_ids)
potential_cross_roads = non_merged_roads.copy()
for road_i in potential_cross_roads :
    non_merged_roads.remove(road_i)
    i = roads_ids.index(road_i)
    node_A = str(roads['NODE_A'][i])
    node_B = str(roads['NODE_B'][i])
    merged = False
    group_number = 0
    # for each group to merge, check if road_i intersects at least two of the group's roads
    while not (merged or group_number == group_number_max):
        if group_number in groups_ids :
            nodes = [node_A, node_B] # list of nodes that have not been identified as connected to the group yet
            for node in group_nodes[group_number] :
                if (str(node) == node_A) and (node_A in nodes) :
                    nodes.remove(node_A)
                elif (str(node) == node_B) and (node_B in nodes) :
                    nodes.remove(node_B)
            if len(nodes) == 0 : # i.e. if both nodes are connected to the group's parallel roads
                merging_roads_groups[group_number].append(road_i)
                Roads_to_groups.loc[i, 'GROUP_ID'] = group_number
                Roads_to_groups.loc[i, 'COLOR'] = group_number % number_of_colors
                Roads_to_groups.loc[i, 'PARALLEL'] = False
                Roads_to_groups.loc[i, 'MERGED'] = True
                Roads_to_groups.loc[i, 'STEP'] = 3
                merged = True
            else :
                group_number = group_number + 1
        else :
            group_number = group_number + 1
    if not merged :
        non_merged_roads.append(road_i)
time_cross_roads_sec = time.time()-start
time_CR = time.gmtime(time_cross_roads_sec)
print("non merged roads without crossroads:", len(non_merged_roads), non_merged_roads)
print("It took", time_CR.tm_hour, 'hour', time_CR.tm_min, 'min', time_CR.tm_sec ,'sec.')

########################## Compute transforms and save #################

# 5. Compute the inverse transform
# Groups = {'GROUP_ID':[], 'ROAD_IDS':[], 'PARALLEL_ROADS_IDS':[], 'NODES_IDS':[], 'MAX_LENGTH':[]}
print('Compute the inverse transform')
max_g = max(groups_ids)
n_street_roads = len(non_merged_roads)
print('number of non merged roads',n_street_roads)
print("groups ids", groups_ids, max_g)
## for identified streets merging roads, save the list of roads to merge
for kth_group in groups_ids:
    Groups[kth_group] = {'ROAD_IDS':[], 'PARALLEL_ROADS_IDS':[], 'NODES_IDS':[], 'MAX_LENGTH':[], 'TUNNEL':[]}
    Groups[kth_group]['NODES_IDS'] = list(set(group_nodes[kth_group]))
    Groups[kth_group]['ROAD_IDS'] = merging_roads_groups[kth_group]
    Groups[kth_group]['PARALLEL_ROADS_IDS'] = parallel_roads[kth_group]
    Groups[kth_group]['MAX_LENGTH'] = groups_max_length[kth_group]
    Groups[kth_group]['TUNNEL'] = tunnel[kth_group]
## for identified streets equal to a road, save their information
for k in range(n_street_roads):
    link_k = non_merged_roads[k]
    # print(Roads_to_groups[Roads_to_groups.ROAD_ID == link_k], link_k)
    road_index = Roads_to_groups[Roads_to_groups.ROAD_ID == link_k].index[0]
    kth_group = max_g + 1 + k
    Roads_to_groups.loc[road_index, 'GROUP_ID'] = kth_group
    Roads_to_groups.loc[road_index, 'PARALLEL'] = True
    Roads_to_groups.loc[road_index, 'MERGED'] = False
    Roads_to_groups.loc[road_index, 'COLOR'] = kth_group % number_of_colors
    node_A = roads['NODE_A'][road_index]
    node_B = roads['NODE_B'][road_index]
    Groups[kth_group] = {'ROAD_IDS': [], 'PARALLEL_ROADS_IDS': [], 'NODES_IDS': [], 'MAX_LENGTH': []}
    Groups[kth_group]['NODES_IDS'] = [node_A, node_B]
    Groups[kth_group]['ROAD_IDS'] = [link_k]
    Groups[kth_group]['PARALLEL_ROADS_IDS'] = [link_k]
    Groups[kth_group]['MAX_LENGTH'] = roads["LENGTH"][road_index]
    Groups[kth_group]['TUNNEL'] = Roads_to_groups.loc[road_index, 'TUNNEL']
# print(Groups)

# 6. Get groups information for nodes
Nodes = {'NODE_ID':[], 'X':[], 'Y':[], 'ROAD_LINKS':[],'GROUPS':[], 'N_GROUPS':[]}
for i, node_i_links in enumerate(nodes_OSM['LINKS']):
    node_i_groups = []
    for link in node_i_links :
        if str(link) in roads_ids :
            link_index = Roads_to_groups[Roads_to_groups.ROAD_ID == str(link)].index[0]
            # print(link_index)
            new_group = Roads_to_groups.loc[link_index, 'GROUP_ID']
            node_i_groups.append(new_group)
                # print("link", link,"link_index", link_index, "group", new_group)
    node_i_groups = list(set(node_i_groups))
    if len(node_i_groups) > 0:
        Nodes['NODE_ID'].append(nodes_OSM['ID'][i])
        Nodes['X'].append(nodes_OSM['X'][i])
        Nodes['Y'].append(nodes_OSM['Y'][i])
        Nodes['ROAD_LINKS'].append(node_i_links)
        Nodes['GROUPS'].append(node_i_groups)
        Nodes['N_GROUPS'].append(len(node_i_groups))
# print(Nodes)

time_total = time.time()
time_t = time.gmtime(time_total)
print("It took", time_t.tm_hour, 'hour', time_t.tm_min, 'min', time_t.tm_sec ,'sec.')

# Save direct and inverse link transform
n_streets = n_street_roads + n_groups
print("Number of streets: ",n_streets)
print("Number of groups: ",n_groups)
pickle.dump(Groups, open(fout_groups, 'wb'))
Roads_to_groups.to_csv(fout_roads, index=False)
pickle.dump(Nodes, open(fout_nodes, 'wb'))

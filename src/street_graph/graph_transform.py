##################################################################################
# Study: Street network generation
# Purpose: Transform the road geo-graph into a street geo-graph
# Note: first generate the street links, then adjust their geometry at the intersection of streets
##################################################################################
import pickle, os, time, sys
import pandas as pd
import numpy as np
import shapely
import skspatial.objects
from skspatial.objects import Line, Point, Points
sys.path.append("../../include")
from lines_geometry import *
from create_poly_around_road import create_1_poly_around_road
from get_indexes import get_indexes
from shapely.geometry import Polygon, LineString
from shapely.ops import nearest_points, transform
import pyproj
start_time = time.time()

# Parameters
angle_threshold_same_line = 30.0 #°
merge_nodes_distance = 10 #meters
width_of_road_buffer = 50 #meters
recover_length_share = 0.8

# Directories
path = "../../"
cache_dir  = path + "temp/"
groups_dir = cache_dir + "street_graph/groups/"
graph_transform_dir = cache_dir + "street_graph/graph_transform/"

# Input and output files
roads_f = cache_dir + "data/load_OSM/roads_links.csv"
groups_f = groups_dir + "groups.dat"
roads_to_groups_f = groups_dir + "roads_to_groups.csv"
nodes_f = groups_dir + 'group_nodes.dat'
street_nodes_f = graph_transform_dir + "street_nodes.csv"
street_links_f = graph_transform_dir + "street_links.csv"
links_transform_f = graph_transform_dir + "links_transform.dat"
nodes_transform_f = graph_transform_dir + "nodes_transform.dat"
assessment_dir = graph_transform_dir + "model_assessment/"
loop_links_f = assessment_dir + "loop_links.csv"
loops_transform_f = assessment_dir + "loops_transform.dat"
links_inverse_transform_f = assessment_dir + "links_inverse_transform.dat"

# Load and format data
if not os.path.exists(graph_transform_dir):
    os.makedirs(graph_transform_dir)
if not os.path.exists(assessment_dir):
    os.makedirs(assessment_dir)
f = open(graph_transform_dir + 'graph_run.txt', 'w')
groups = pickle.load(open(groups_f, 'rb'))
groups_df = pd.DataFrame.from_dict(groups, orient='index')
groups_df = groups_df.reset_index(names="GROUP_ID")
road_nodes = pickle.load(open(nodes_f,'rb'))
Road_nodes = pd.DataFrame(road_nodes)
Road_nodes = Road_nodes.astype({'NODE_ID':str})
nodes_list = list(road_nodes['NODE_ID'])
n_road_nodes = len(nodes_list)
roads = pd.read_csv(roads_f)
roads['GEOMETRY'] = roads['GEOMETRY'].apply(lambda x: shapely.wkt.loads(x))
print("Data loaded...")

# Initialize output dataframes
geo_transformer = pyproj.Transformer.from_crs('EPSG:2154', 'EPSG:3857', always_xy=True)
Street_links = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'TUNNEL'])
Groups_nodes = pd.DataFrame(columns=['GROUP_NODE', 'X', 'Y','GROUP','GROUP_SLOPE','ROAD_NODES'])
links_transform_df = pd.read_csv(roads_to_groups_f)
links_transform_df['GEOMETRY'] = links_transform_df['GEOMETRY'].apply(lambda x: shapely.wkt.loads(x))
n_roads = links_transform_df.shape[0]
links_transform_df['STREET_IDS']=[[] for r in range(n_roads)]
links_transform_df['SHARES']=[[] for s in range(n_roads)] # links_transform_df columns=['ROAD_ID', 'GROUP_ID', 'MAIN_ROAD', 'MERGED', 'COLOR', 'GEOMETRY']) + add 'STREET_IDS'and 'SHARES'
links_transform = links_transform_df.to_dict(orient='list')
road_ids = [str(road) for road in roads['ID']]
street_to_roads = {}
nodes_transform = {}
for road_node in nodes_list :
    nodes_transform[str(road_node)] = {'GROUP_NODES':[], 'STREET_NODE':0}
print("Outputs initialized...")

# Step 1 : Get roads and links within a group
line="Step 1 : Get roads and links within a group\n"
print(line)
f.write(line)
nodes_cnt = 1
links_cnt = 1
for i, group_roads in enumerate(groups_df['MAIN_ROAD_IDS']) :
    if i % 1000 == 0 :
        time_street_links = time.time() - start_time
        time_gm = time.gmtime(time_street_links)
        line = str(i) + "streets took"+ str(time_gm.tm_hour) + 'hour' + str(time_gm.tm_min) + 'min' + str(time_gm.tm_sec) + 'sec.'
        print(line)
        f.write(line + '\n')
    tunnel = groups_df.loc[i, "TUNNEL"]
    group_i = groups_df.loc[i, "GROUP_ID"]
    if len(group_roads) > 1 :
        # a) Among main roads of the group, identify those on the same line
        lines = {'ROAD_IDS':[], 'POINT_A':[],'POINT_B':[]}
        for road_id in group_roads :
            road_idx = road_ids.index(road_id)
            road_xa = roads.loc[road_idx, 'XA']
            road_ya = roads.loc[road_idx, 'YA']
            road_xb = roads.loc[road_idx, 'XB']
            road_yb = roads.loc[road_idx, 'YB']
            pointA = [road_xa, road_ya]
            pointB = [road_xb, road_yb]
            n_tests = len(lines['ROAD_IDS'])
            # if it is the first road of the group, just add the associated line
            if n_tests == 0 :
                lines['ROAD_IDS'].append([road_id])
                lines['POINT_A'].append(pointA)
                lines['POINT_B'].append(pointB)
            # else, check if the road belongs to the same line as another road of the group
            else :
                existing_line = False
                k = 0
                while not (existing_line or k == n_tests):
                    pointC = lines['POINT_A'][k]
                    pointD = lines['POINT_B'][k]
                    existing_line = segments_on_the_same_line(pointA, pointB, pointC, pointD, angle_threshold_same_line)
                    if existing_line : # if it does, add its ID to the line
                        lines['ROAD_IDS'][k].append(road_id)
                    else :
                        k = k + 1
                # else, create a new line
                if not existing_line :
                    lines['ROAD_IDS'].append([road_id])
                    lines['POINT_A'].append(pointA)
                    lines['POINT_B'].append(pointB)
        # b) Define the line of the group
        n_lines = len(lines['ROAD_IDS'])
        if n_lines > 1 :
            node_to_project = lines['POINT_A'][0]
            s, y_inter = get_street_line(node_to_project, lines['POINT_A'], lines['POINT_B'])
            street_line = Line.from_slope(slope=s,y_intercept=y_inter)
        else :
            s = slope([lines['POINT_A'][0][0], lines['POINT_A'][0][1], lines['POINT_B'][0][0], lines['POINT_B'][0][1]])
            y_inter = lines['POINT_A'][0][1] - (s * lines['POINT_A'][0][0])
            street_line = Line.from_slope(slope=s,y_intercept=y_inter)
        # c) Project all nodes on the street line and merge them if they are too close
        group_nodes = groups[group_i]['NODES_IDS']
        projected_nodes = {'GROUP_NODES':[], 'COORDINATES':[], 'ROAD_NODES':[], 'X':[], 'Y':[]}
        ## first, create the new street nodes within the group
        for node_N in group_nodes :
            node_N = str(node_N)
            node_index = nodes_list.index(node_N)
            x = road_nodes['X'][node_index]
            y = road_nodes['Y'][node_index]
            [xp, yp] = street_line.project_point([x,y])
            pointP = skspatial.objects.Point([xp, yp])
            existing_node = False
            ### check if it merges with previously created group_nodes
            n_test = len(projected_nodes['GROUP_NODES'])
            k = 0
            while (not existing_node) and (k < n_test):
                coordinates = projected_nodes['COORDINATES'][k][0]
                distance = pointP.distance_point(coordinates)
                if distance < merge_nodes_distance : # if there are less than 10 meters between two nodes, they are the same
                    existing_node = True
                    projected_nodes['ROAD_NODES'][k].append(node_N)
                    projected_nodes['COORDINATES'][k].append([xp, yp])
                else :
                    k = k+1
            if not existing_node : ### if not, create a new group_node
                projected_nodes['GROUP_NODES'].append(nodes_cnt)
                projected_nodes['COORDINATES'].append([[xp,yp]])
                projected_nodes['ROAD_NODES'].append([node_N])
                nodes_cnt += 1
        ## second, if several nodes are merged, get their centroid coordinates
        for node_coordinates_list in projected_nodes['COORDINATES']:
            if len(node_coordinates_list) > 1:
                projected_points = Points(node_coordinates_list)
                points_centered, centroid = projected_points.mean_center(return_centroid=True)
                projected_nodes['X'].append(centroid[0])
                projected_nodes['Y'].append(centroid[1])
            else:
                projected_nodes['X'].append(node_coordinates_list[0][0])
                projected_nodes['Y'].append(node_coordinates_list[0][1])
        ## third, add close junction nodes to the road-nodes lists
        junction_nodes = groups[group_i]["JUNCTION_NODES"]
        street_nodes_list = list(zip(projected_nodes['X'], projected_nodes['Y']))
        street_nodes = shapely.MultiPoint(street_nodes_list)
        for node_N in junction_nodes :
            node_N = str(node_N)
            node_index = nodes_list.index(node_N)
            x = road_nodes['X'][node_index]
            y = road_nodes['Y'][node_index]
            pointP = shapely.Point([x, y])
            matching_point_coord = nearest_points(pointP, street_nodes)[1]
            matching_point = street_nodes_list.index((matching_point_coord.x, matching_point_coord.y))
            if node_N not in projected_nodes['ROAD_NODES'][matching_point] :
                projected_nodes['ROAD_NODES'][matching_point].append(node_N)
        ## fourth, save information in Groups_nodes
        Group_i_nodes = pd.DataFrame(columns=['GROUP_NODE', 'X', 'Y','GROUP','GROUP_SLOPE','ROAD_NODES'])
        for node_idx, node_id in enumerate(projected_nodes['GROUP_NODES']):
            new_node = pd.DataFrame(columns=['GROUP_NODE', 'X', 'Y','GROUP','GROUP_SLOPE','ROAD_NODES'])
            new_node.loc[0, 'GROUP_NODE'] = node_id
            new_node.loc[0, 'X'] = projected_nodes['X'][node_idx]
            new_node.loc[0, 'Y'] = projected_nodes['Y'][node_idx]
            new_node.loc[0, 'GROUP'] = group_i #groups[group_i]['GROUP_ID']
            new_node.loc[0, 'GROUP_SLOPE'] = s
            matching_road_nodes = projected_nodes['ROAD_NODES'][node_idx]
            new_node.loc[0, 'ROAD_NODES'] = matching_road_nodes
            Group_i_nodes = pd.concat([Group_i_nodes, new_node], ignore_index=True)
            for road_node in matching_road_nodes :
                nodes_transform[str(road_node)]['GROUP_NODES'].append(node_id)
                nodes_transform[str(road_node)]['STREET_NODE'] = node_id
        ## order nodes to prepare the generation of links
        Group_i_nodes = Group_i_nodes.sort_values(by=['X']).reset_index(drop=True)
        n_nodes = Group_i_nodes.shape[0]
        if (n_nodes > 1) and (Group_i_nodes.loc[0, 'X'] == Group_i_nodes.loc[n_nodes - 1, 'X']):
            Group_i_nodes = Group_i_nodes.sort_values(by=['Y']).reset_index(drop=True)
        Groups_nodes = pd.concat([Groups_nodes, Group_i_nodes], ignore_index=True)
        # d) Generate street links and association between road_links and street_links
        Group_links = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'AXIAL_LINE','LENGTH', 'TUNNEL'])
        ## Create street links within the group
        if n_nodes > 1 :
            for j in range(n_nodes-1):
                new_link = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'AXIAL_LINE','LENGTH','TUNNEL'])
                new_link.loc[0,'STREET_ID'] = links_cnt
                new_link.loc[0, 'NODE_A'] = Group_i_nodes.loc[j,'GROUP_NODE']
                new_link.loc[0, 'NODE_B'] = Group_i_nodes.loc[j+1,'GROUP_NODE']
                new_link.loc[0, 'TUNNEL'] = tunnel
                xa = Group_i_nodes.loc[j, 'X']
                ya = Group_i_nodes.loc[j, 'Y']
                xb = Group_i_nodes.loc[j+1, 'X']
                yb = Group_i_nodes.loc[j+1, 'Y']
                detailed_geometry = LineString([[xa,ya],[xb,yb]])
                new_link.loc[0, 'AXIAL_LINE'] = detailed_geometry
                line_transformed = transform(geo_transformer.transform, detailed_geometry)
                new_link.loc[0, 'LENGTH'] = line_transformed.length
                links_cnt = links_cnt + 1
                Group_links = pd.concat([Group_links, new_link], ignore_index=True)
        else :
            new_link = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'AXIAL_LINE','LENGTH', 'TUNNEL'])
            new_link.loc[0, 'STREET_ID'] = links_cnt
            new_link.loc[0, 'NODE_A'] = Group_i_nodes.loc[0, 'GROUP_NODE']
            new_link.loc[0, 'NODE_B'] = Group_i_nodes.loc[0, 'GROUP_NODE']
            new_link.loc[0, 'TUNNEL'] = tunnel
            xa = Group_i_nodes.loc[0, 'X']
            ya = Group_i_nodes.loc[0, 'Y']
            new_link.loc[0, 'AXIAL_LINE'] = LineString([[xa, ya], [xa, ya]])
            new_link.loc[0, 'LENGTH'] = 0
            links_cnt = links_cnt + 1
            Group_links = pd.concat([Group_links, new_link], ignore_index=True)
        ## Associate main road_links with street_links within the group
        for road_r in group_roads :
            r = road_ids.index(road_r)
            xa = roads.loc[r,'XA']
            ya = roads.loc[r,'YA']
            xb = roads.loc[r,'XB']
            yb = roads.loc[r,'YB']
            road_length = roads.loc[r,'LENGTH']
            road_poly = create_1_poly_around_road(xa,ya,xb,yb,width_of_road_buffer)
            street_ids = []
            shares = []
            for s, street_s in enumerate(Group_links['AXIAL_LINE']):
                if road_poly.intersects(street_s) :
                    intersection = road_poly.intersection(street_s).length
                    street_length = street_s.length
                    if (intersection > recover_length_share * street_length) or (intersection > recover_length_share * road_length) :
                        identified_street = Group_links.loc[s,'STREET_ID']
                        street_ids.append(identified_street)
                        shares.append(street_length) # the initial road can contribute to several streets if it is long
                        if identified_street in street_to_roads.keys() :
                            street_to_roads[identified_street].append(road_r)
                        else :
                            street_to_roads[identified_street] = [road_r]
            if len(street_ids) > 0 :
                total_streets_length = np.sum(shares)
                links_transform['STREET_IDS'][r].extend(street_ids)
                links_transform['SHARES'][r].extend([share / total_streets_length for share in shares])
            else : #  if no street identified, just identify close links
                pointP = shapely.Point([roads.loc[r, 'XA'], roads.loc[r, 'YA']])
                closest_link_index = closest_line_from_point(pointP, list(Group_links['AXIAL_LINE']))
                closest_link = Group_links.loc[closest_link_index, 'STREET_ID']
                if closest_link in street_to_roads.keys() :
                    street_to_roads[closest_link].append(road_r)
                else:
                    street_to_roads[closest_link] = [road_r]
                links_transform['STREET_IDS'][r].append(closest_link)
                links_transform['SHARES'][r].append(1)
        ## Associate junction links with street_links within the group
        junctions = groups[group_i]["JUNCTION_IDS"]
        for road_j in junctions : ### look for the nearest street link
            r = road_ids.index(road_j)
            pointC = shapely.Point([roads.loc[r,'XA'], roads.loc[r,'YA']])
            closest_link_index = closest_line_from_point(pointC, list(Group_links['AXIAL_LINE']))
            closest_link = Group_links.loc[closest_link_index,'STREET_ID']
            if closest_link in street_to_roads.keys() :
                street_to_roads[closest_link].append(road_j)
            else :
                street_to_roads[closest_link] = [road_j]
            links_transform['STREET_IDS'][r].append(closest_link)
            links_transform['SHARES'][r].append(1)
        Group_links = Group_links[['STREET_ID', 'NODE_A', 'NODE_B', 'TUNNEL', 'AXIAL_LINE', 'LENGTH']]
        Street_links = pd.concat([Street_links, Group_links], ignore_index=True)
    else : # non-merged roads
        # Node id generation and save node transform
        Group_i_nodes = pd.DataFrame(columns=['GROUP_NODE', 'X', 'Y', 'GROUP', 'GROUP_SLOPE', 'ROAD_NODES'])
        [road_node_A,road_node_B] = groups[group_i]['NODES_IDS']
        node_A_index = nodes_list.index(str(road_node_A))
        node_B_index = nodes_list.index(str(road_node_B))
        xa = road_nodes['X'][node_A_index]
        ya = road_nodes['Y'][node_A_index]
        xb = road_nodes['X'][node_B_index]
        yb = road_nodes['Y'][node_B_index]
        s = slope([xa,ya,xb,yb])
        ## add node_A to Groups_nodes
        Group_i_nodes.loc[0, 'GROUP_NODE'] = nodes_cnt
        Group_i_nodes.loc[0, 'ROAD_NODES'] = [str(road_node_A)]
        Group_i_nodes.loc[0, 'X'] = xa
        Group_i_nodes.loc[0, 'Y'] = ya
        Group_i_nodes.loc[0, 'GROUP'] = group_i
        Group_i_nodes.loc[0, 'GROUP_SLOPE'] = s
        ## add node_B to Groups_nodes
        Group_i_nodes.loc[1, 'GROUP_NODE'] = nodes_cnt + 1
        Group_i_nodes.loc[1, 'ROAD_NODES'] = [str(road_node_B)]
        Group_i_nodes.loc[1, 'X'] = xb
        Group_i_nodes.loc[1, 'Y'] = yb
        Group_i_nodes.loc[1, 'GROUP'] = group_i
        Group_i_nodes.loc[1, 'GROUP_SLOPE'] = s
        Groups_nodes = pd.concat([Groups_nodes, Group_i_nodes], ignore_index=True)
        ## nodes tranform
        nodes_transform[str(road_node_A)]['GROUP_NODES'].append(nodes_cnt)
        nodes_transform[str(road_node_A)]['STREET_NODE'] = nodes_cnt
        nodes_transform[str(road_node_B)]['GROUP_NODES'].append(nodes_cnt+1)
        nodes_transform[str(road_node_B)]['STREET_NODE'] = nodes_cnt+1
        ## link transform
        road_r = group_roads[0]
        r = road_ids.index(road_r)
        links_transform['STREET_IDS'][r] = [links_cnt]
        links_transform['SHARES'][r] = [1]
        ## link in Street_links
        new_link = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'TUNNEL', 'AXIAL_LINE', 'LENGTH'])
        detailed_geometry = links_transform_df.loc[r,'GEOMETRY']
        new_link.loc[0, 'STREET_ID'] = links_cnt
        new_link.loc[0, 'NODE_A'] = nodes_cnt
        new_link.loc[0, 'NODE_B'] = nodes_cnt + 1
        new_link.loc[0, 'TUNNEL'] = tunnel
        new_link.loc[0, 'AXIAL_LINE'] = detailed_geometry
        line_transformed = transform(geo_transformer.transform, detailed_geometry)
        new_link.loc[0, 'LENGTH'] = line_transformed.length
        Street_links = pd.concat([Street_links, new_link], ignore_index=True)
        if links_cnt in street_to_roads.keys():
            street_to_roads[links_cnt].append(road_r)
        else:
            street_to_roads[links_cnt] = [road_r]
        links_cnt = links_cnt + 1
        nodes_cnt = nodes_cnt + 2

# Step 2 : Calculate coordinates of the groups summits at the intersection of streets using Groups_nodes
line="Step 2 : Calculate coordinates of the groups summits at the intersection of streets"
print(line)
f.write(line + '\n')
start_inter = time.time()
Intersection_nodes = Road_nodes[(Road_nodes['N_GROUPS'] > 1)]
Groups_nodes = Groups_nodes.rename(columns={'GROUP_NODE':'NODE_ID'})
group_nodes_ids = list(Groups_nodes['NODE_ID'])
intersection_road_nodes = list(Intersection_nodes['NODE_ID'])
line = "number of intersections: " + str(len(intersection_road_nodes))
print(line)
f.write(line + '\n')
nodes_idx_to_exclude = []
while len(intersection_road_nodes) > 0 :
    if len(intersection_road_nodes) % 1000 == 0 :
        time_inter_sec = time.time() - start_inter
        time_tot = time.gmtime(time_inter_sec)
        line = str(len(intersection_road_nodes)) + "intersections took"+ str(time_tot.tm_hour) + 'hour' + str(time_tot.tm_min) + 'min' + str(
            time_tot.tm_sec) + 'sec.'
        print(line)
        f.write(line + '\n')
    road_node_i = intersection_road_nodes[0]
    road_node_index = nodes_list.index(road_node_i)
    node_i_groups = road_nodes['GROUPS'][road_node_index]
    # if the node is an intersection of streets, calculate its coordinates and replace group_node by street_node in transforms
    if len(node_i_groups) > 1 :
        # 1) Calculate the coordinates of the intersection of the groups/streets
        ## get information on group nodes corresponding to the road_node
        Groups_info = pd.DataFrame(columns = Groups_nodes.columns)
        cnt = 0
        road_nodes_test = [road_node_i]
        merging_road_nodes = []
        potential_group_nodes = []
        while any(group_node not in list(Groups_info['NODE_ID']) for group_node in potential_group_nodes) or cnt == 0: 
            merging_road_nodes = merging_road_nodes + road_nodes_test
            road_nodes_test_next_loop = []
            for road_node_j in road_nodes_test:  # for each road node Nj to test
                matching_group_nodes_j = nodes_transform[str(road_node_j)]['GROUP_NODES']
                for group_node_J in matching_group_nodes_j : # for each of its groups,
                    # if it was not checked previously, look for the group node matching Nj
                    J_index = group_nodes_ids.index(group_node_J)
                    new_row = Groups_nodes[Groups_nodes.index == J_index].reset_index(drop=True)
                    if group_node_J not in list(Groups_info['NODE_ID']):
                        Groups_info = pd.concat([Groups_info, new_row], ignore_index=True)
                    # and save the groups of other road nodes matching Nj'
                    road_nodes_projected_on_J = new_row.loc[0, 'ROAD_NODES']
                    road_nodes_test_next_loop = list(set(road_nodes_test_next_loop + [node for node in road_nodes_projected_on_J if node not in merging_road_nodes]))
                road_nodes_test = road_nodes_test_next_loop
            for road_node_k in road_nodes_test_next_loop:
                potential_group_nodes.extend(nodes_transform[str(road_node_k)]['GROUP_NODES'])
            cnt += 1
        merging_road_nodes = merging_road_nodes + road_nodes_test
        for checked_road_node in merging_road_nodes:  # remove identified intersection road nodes
            if str(checked_road_node) in intersection_road_nodes: ### str(checked_road_node) not necessarily in intersection_road_nodes
                intersection_road_nodes.remove(str(checked_road_node))
        ## get the intersection of pairs of roads
        n_groups = len(Groups_info)
        intersection_x = []
        intersection_y = []
        intersection_coord = []
        for i in range(n_groups):
            x = Groups_info.loc[i, 'X']
            y = Groups_info.loc[i, 'Y']
            if (x not in intersection_x) or (y not in intersection_y):
                intersection_coord.append([x, y])
                intersection_x.append(x)
                intersection_y.append(y)
        ## get the centroid of the intersections
        if len(intersection_coord) == 1 :
            centroid = intersection_coord[0]
        else :
            points_centered, centroid = Points(intersection_coord).mean_center(return_centroid=True)
        xi = centroid[0]
        yi = centroid[1]
        # 2) Modify the nodes information (id and coordinates) where they appear
        street_node = Groups_info.loc[0,'NODE_ID'] # select the first merging group_node: we will take its id as street id
        road_nodes_in_group = []
        for k, group_node in enumerate(Groups_info['NODE_ID']) :
            ## first, in the global list of nodes
            group_node_index = group_nodes_ids.index(group_node)
            if k == 0 : # then, the column 'GROUP_NODE' will be labelled 'NODE_ID', refering to the street node id
                Groups_nodes.loc[group_node_index,'NODE_ID'] = street_node
                Groups_nodes.loc[group_node_index, 'X'] = xi
                Groups_nodes.loc[group_node_index, 'Y'] = yi
            else :
                nodes_idx_to_exclude.append(group_node_index)
            ## add id to nodes_transform
            road_nodes_in_group = list(set(road_nodes_in_group + Groups_nodes.loc[group_node_index,'ROAD_NODES']))
            ## modify the node is where it appears in Street_links
            nodes_A_to_modify = get_indexes(list(Street_links['NODE_A']),group_node)
            if len(nodes_A_to_modify) > 0 :
                for node in nodes_A_to_modify :
                    Street_links.loc[node,'NODE_A'] = street_node
            nodes_B_to_modify = get_indexes(list(Street_links['NODE_B']), group_node)
            if len(nodes_B_to_modify) > 0 :
                for node in nodes_B_to_modify:
                    Street_links.loc[node, 'NODE_B'] = street_node
        for road_node_i in road_nodes_in_group :
            nodes_transform[str(road_node_i)]['STREET_NODE'] = street_node # remark : this is the only case where the street_node ID is different from the group_node_ID
    else : # if it is not an intersection
        line = "not an intersection: node "+ str(road_node_i)+ "- road links: "+ str(road_nodes['ROAD_LINKS'][road_node_index])
        print(line)
        f.write(line)
Groups_nodes = Groups_nodes[~Groups_nodes.index.isin(nodes_idx_to_exclude)]
time_inter_sec = time.time() - start_inter
time_tot = time.gmtime(time_inter_sec)
line = str(len(intersection_road_nodes)) + "intersections took"+ str(time_tot.tm_hour) + 'hour' + str(time_tot.tm_min) + 'min' + str(
            time_tot.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')

# Step 3) Generate street links geometry
line = "Step 3) Generate street links geometry"
print(line)
f.write(line + '\n')
start_geom = time.time()
Street_nodes = Groups_nodes[['NODE_ID', 'X', 'Y']]
n_links = Street_links.shape[0]
links_to_exclude = {}
for i, street_link in enumerate(Street_links['STREET_ID']):
    # for each link, get its nodes ids
    node_A_id = Street_links.loc[i,'NODE_A']
    node_B_id = Street_links.loc[i,'NODE_B']
    if node_A_id == node_B_id :
        links_to_exclude[street_link] = Street_links.loc[i,'NODE_A']
    # then, their geographical position
    node_A = Street_nodes[Street_nodes.NODE_ID.isin([int(node_A_id),str(node_A_id)])].reset_index(drop=True)
    node_B = Street_nodes[Street_nodes.NODE_ID.isin([int(node_B_id),str(node_B_id)])].reset_index(drop=True)
    xa = node_A.loc[0,'X']
    ya = node_A.loc[0, 'Y']
    xb = node_B.loc[0, 'X']
    yb = node_B.loc[0, 'Y']
    # finally, add the link geometry
    Street_links.loc[i, 'GEOMETRY'] = LineString([[xa,ya],[xb,yb]])
    Street_links.loc[i, 'XA'] = xa
    Street_links.loc[i, 'YA'] = ya
    Street_links.loc[i, 'XB'] = xb
    Street_links.loc[i, 'YB'] = yb
street_transformed = list(street_to_roads.keys())
Street_links = Street_links[Street_links.STREET_ID.isin(street_transformed)] # filter streets without associated roads
line = str(len(links_to_exclude.keys())) + "links to exclude"
print(line)
f.write(line + '\n')
time_inter_sec = time.time() - start_geom
time_tot = time.gmtime(time_inter_sec)
print("Geometry took", time_tot.tm_hour, 'hour', time_tot.tm_min, 'min', time_tot.tm_sec, 'sec.')

### ici, récupérer les road_links et les transférer à un autre street_link avec nodes en commun.
line="Step 3b) Filter loop-links"
print(line)
f.write(line + '\n')
Loop_links = Street_links[Street_links.STREET_ID.isin(list(links_to_exclude.keys()))]
loops_to_roads={}
for link in Loop_links['STREET_ID']:
    if int(link) in street_to_roads.keys():
        loops_to_roads[int(link)] = street_to_roads[int(link)]
Street_links = Street_links[~Street_links.STREET_ID.isin(list(links_to_exclude.keys()))]
for street_link, street_node in links_to_exclude.items() :
    links_connected_to_node = Street_links[(Street_links.NODE_A == street_node) | (Street_links.NODE_B == street_node)]
    links_connected_to_node_filtered = links_connected_to_node[~((links_connected_to_node.NODE_A == street_node)&(Street_links.NODE_B == street_node))]
    filtered_links = list(links_connected_to_node_filtered['STREET_ID'])
    matching_link = filtered_links[0]
    # replace the link id in links_transform
    if street_link in street_to_roads.keys() :
        roads_to_move = street_to_roads[street_link] ### if this is to save, change roads_to_move key
        for road_r in roads_to_move :
            r = road_ids.index(road_r)
            projection_on_streets = links_transform['STREET_IDS'][r]
            index_to_replace = projection_on_streets.index(int(street_link))
            links_transform['STREET_IDS'][r][index_to_replace] = matching_link
nodes_transform_df = pd.DataFrame.from_dict(nodes_transform,orient='index').reset_index(drop=False)
time_tot_sec = time.time() - start_time
time_tot = time.gmtime(time_tot_sec)
line = "Total time:"+ str(time_tot.tm_mday) + 'days'+ str(time_tot.tm_hour) + 'hour' + str(time_tot.tm_min) + 'min' + str(
            time_tot.tm_sec) + 'sec.'
print(line)
f.write(line + '\n')
print("Number of streets:", Street_links.shape[0])
print("Number of nodes:", Street_nodes.shape[0])

# Save files as csv
Street_nodes.to_csv(street_nodes_f)
Street_links.to_csv(street_links_f)
pickle.dump(links_transform, open(links_transform_f,'wb'))
pickle.dump(nodes_transform_df, open(nodes_transform_f,'wb'))
Loop_links.to_csv(loop_links_f)
pickle.dump(loops_to_roads, open(loops_transform_f,'wb'))
pickle.dump(street_to_roads, open(links_inverse_transform_f,'wb'))
f.close()

##################################################################################
# Study: Street network generation
# Purpose: Transform the road geo-graph into a street geo-graph
# Author: Marjolaine Lannes
# Creation date: April 25, 2023
# Note: first generate the street links, then adjust their geometry at streets intersections
##################################################################################
## AJOUTER UNE CONDITION : si les deux noeuds d'une arête sont connectés à  un même groupe :
# l'arrête s'ajoute au groupe. +++ tunnel
## Penser à avoir les road_ids sous str pas int : c'est le cas dans Roads_to_groups
import pickle
import pandas as pd
import numpy as np
from skspatial.objects import Line, Point, Points
from include.lines_geometry import *
from include.create_poly_around_road import create_1_poly_around_road
from include.get_indexes import *
from shapely.geometry import Polygon, LineString
# project_point_on_line / distance_to_line / lines_intersection / segments_on_the_same_line / get_street_line

# Input and output files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
#groups_dir = path + "temp/street_graph/groups/"
groups_dir = path + "temp/street_graph/groups_30_08/"
groups_f = groups_dir + "groups.dat"
roads_to_groups_f = groups_dir + "roads_to_groups.csv"
nodes_f = groups_dir + 'group_nodes.dat'
# dir_roads = path + "temp/load_OSM/subdomains_roads_with_buffer_100m_lamb93/"
dir_roads = path + "temp/load_OSM/subdomains_roads_lamb93/"
outdir = path + "temp/street_graph/graph_transform/"
street_nodes_f = outdir + "street_nodes.csv"
street_links_f = outdir + "street_links.csv"
street_links_dat_f = outdir + "street_links.dat"
links_transform_f = outdir + "links_transform.dat"
nodes_transform_f = outdir + "nodes_transform.dat"

# Parameters
distance_merge_points = 15 # meters

# Cell data
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
grid = pickle.load(open(gridfile,'rb'))
poly_num = 4672
x0 = grid['we_id'][poly_num]
y0 = grid['sn_id'][poly_num]
prefix='we_id_from_' + str(x0) + '_sn_id_from_' + str(y0)
fname = 'sub_' + prefix + '.dat'

# Load data
groups = pickle.load(open(groups_f, 'rb'))
# groups[i] = {'ROAD_IDS':[], 'PARALLEL_ROADS_IDS':[], 'NODES_IDS':[] ,'LONGEST_ROAD':[], 'MAX_LENGTH':[]}
print(groups[45])
road_nodes = pickle.load(open(nodes_f,'rb'))
nodes_list = [str(node) for node in road_nodes['NODE_ID']]
road_nodes['NODE_ID'] = nodes_list
Road_nodes = pd.DataFrame(road_nodes)
roads = pickle.load(open(dir_roads + fname, 'rb')) ### CELL -> GLOBAL
roads['ROAD_ID'] = [str(road) for road in roads['ROAD_ID']]
print("Data loaded...")

# Initialize output dataframes
group_node_col = ['GROUP_NODE', 'X', 'Y','GROUP','GROUP_SLOPE','ROAD_NODES']
Street_links = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY_ATTRIBUTES'])
Groups_nodes = pd.DataFrame(columns=group_node_col)
links_transform_df = pd.read_csv(roads_to_groups_f)
# columns=['ROAD_ID', 'GROUP_ID', 'PARALLEL', 'MERGED', 'COLOR', 'GEOMETRY']) + add 'STREET_IDS'and 'SHARES'
links_transform = links_transform_df.to_dict('list') # set_index('ROAD_ID')
n_roads = links_transform_df.shape[0]
links_transform['STREET_IDS']=[[] for r in range(n_roads)]
links_transform['SHARES']=[[] for s in range(n_roads)]
print(len(links_transform['SHARES']))
#print(links_transform) r = road_ids.index(road_r)
group_nodes_transform = {}
for node in nodes_list :
    group_nodes_transform[node] = {'GROUP_NODES':{}, 'STREET_NODE':None}
nodes_transform = pd.DataFrame(columns=['ROAD_NODE', 'STREET_NODE'])
nodes_transform['ROAD_NODE'] = nodes_list

# Step 1 : Get roads and links within a group
road_ids = roads['ROAD_ID']
nodes_cnt = 1
links_cnt = 1
for group_i, group_info in groups.items() :
    print(group_i)
    group_roads = group_info['PARALLEL_ROADS_IDS']
    if len(group_roads) > 1 :
        # a) Identify parallel roads within the group
        lines = {'ROAD_IDS':[], 'POINT_A':[],'POINT_B':[]}
        for road_id in group_roads :
            i = road_ids.index(road_id)
            road_xa = roads['XA'][i]
            road_ya = roads['YA'][i]
            road_xb = roads['XB'][i]
            road_yb = roads['YB'][i]
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
                    existing_line = segments_on_the_same_line(pointA, pointB, pointC, pointD)
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
        road_nodes_in_group = group_info['NODES_IDS']
        projected_nodes = {'GROUP_NODES':[], 'COORDINATES':[], 'ROAD_NODES':[],'X':[],'Y':[]}
        ## first, create the new street nodes within the group
        for node_Nr in road_nodes_in_group :
            node_Nr = str(node_Nr)
            node_index = nodes_list.index(node_Nr)
            x = road_nodes['X'][node_index]
            y = road_nodes['Y'][node_index]
            [xp, yp] = street_line.project_point([x,y])
            pointP = Point([xp, yp])
            existing_node = False
            ### check if it merges with previously created group_nodes
            n_test = len(projected_nodes['GROUP_NODES'])
            k = 0
            while (not existing_node) and (k < n_test):
                coordinates = projected_nodes['COORDINATES'][k][0]
                distance = pointP.distance_point(coordinates)
                if distance < distance_merge_points : # if there are less than d meters between two nodes, they are the same
                    existing_node = True
                    projected_nodes['ROAD_NODES'][k].append(node_Nr)
                    projected_nodes['COORDINATES'][k].append([xp, yp])
                    group_nodes_transform[node_Nr]['GROUP_NODES'][group_i] = projected_nodes['GROUP_NODES'][k]
                    group_nodes_transform[node_Nr]['STREET_NODE'] = projected_nodes['GROUP_NODES'][k]
                else :
                    k = k+1
            ### if not, create a new group_node
            if not existing_node :
                projected_nodes['GROUP_NODES'].append(nodes_cnt)
                projected_nodes['COORDINATES'].append([[xp,yp]])
                projected_nodes['ROAD_NODES'].append([node_Nr])
                group_nodes_transform[node_Nr]['GROUP_NODES'][group_i] = nodes_cnt
                group_nodes_transform[node_Nr]['STREET_NODE'] = nodes_cnt
                nodes_cnt += 1
        ## second, if several nodes are merged, get their centroid coordinates
        for node_coordinates_list in projected_nodes['COORDINATES']:
            if len(node_coordinates_list)>1 :
                projected_points = Points(node_coordinates_list)
                points_centered, centroid = projected_points.mean_center(return_centroid=True)
                projected_nodes['X'].append(centroid[0])
                projected_nodes['Y'].append(centroid[1])
            else :
                projected_nodes['X'].append(node_coordinates_list[0][0])
                projected_nodes['Y'].append(node_coordinates_list[0][1])
        ## third, save information in Groups_nodes
        Group_i_nodes = pd.DataFrame(columns=group_node_col)
        for i, node_id in enumerate(projected_nodes['GROUP_NODES']):
            new_node = {}
            new_node['GROUP_NODE'] = node_id
            new_node['X'] = projected_nodes['X'][i]
            new_node['Y'] = projected_nodes['Y'][i]
            new_node['GROUP'] = int(group_i)
            new_node['GROUP_SLOPE'] = s
            new_node['ROAD_NODES'] = [projected_nodes['ROAD_NODES'][i]]
            new_node = pd.DataFrame.from_dict(new_node)
            Group_i_nodes = pd.concat([Group_i_nodes, new_node], ignore_index=True)
        Group_i_nodes = Group_i_nodes.sort_values(by=['X']).reset_index(drop=True)
        n_nodes = Group_i_nodes.shape[0]
        if (n_nodes > 1) and (Group_i_nodes.loc[0,'X'] == Group_i_nodes.loc[n_nodes-1,'X']):
            Group_i_nodes = Group_i_nodes.sort_values(by=['Y']).reset_index(drop=True)
        Groups_nodes = pd.concat([Groups_nodes, Group_i_nodes], ignore_index=True)
        # d) Generate street links and association between road_links and street_links
        Group_links = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY'])
        ## Create street links within the group
        for j in range(n_nodes-1):
            new_link = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY'])
            new_link.loc[0,'STREET_ID'] = links_cnt
            new_link.loc[0, 'NODE_A'] = Group_i_nodes.loc[j,'GROUP_NODE']
            new_link.loc[0, 'NODE_B'] = Group_i_nodes.loc[j+1,'GROUP_NODE']
            xa = Group_i_nodes.loc[j, 'X']
            ya = Group_i_nodes.loc[j, 'Y']
            xb = Group_i_nodes.loc[j+1, 'X']
            yb = Group_i_nodes.loc[j+1, 'Y']
            new_link.loc[0, 'GEOMETRY'] = LineString([[xa,ya],[xb,yb]])
            links_cnt = links_cnt + 1
            Group_links = pd.concat([Group_links, new_link], ignore_index=True)
        ## Associate road_links with street_links within the group
        for road_r in group_roads :
            r = road_ids.index(road_r)
            xa = roads['XA'][r]
            ya = roads['YA'][r]
            xb = roads['XB'][r]
            yb = roads['YB'][r]
            road_length = roads['LENGTH'][r]
            width = 50
            road_poly = create_1_poly_around_road(xa,ya,xb,yb,width)
            street_ids = []
            shares = []
            for s, street_s in enumerate(Group_links['GEOMETRY']):
                if road_poly.intersects(street_s) :
                    intersection = road_poly.intersection(street_s).length
                    street_length = street_s.length
                    if intersection > 0.8 * street_length :
                        identified_street = Group_links.loc[s,'STREET_ID']
                        street_ids.append(identified_street)
                        shares.append(street_length/road_length) # the initial road can contribute to several streets if it is long
            #print(street_ids)
            links_transform['STREET_IDS'][r] = street_ids
            links_transform['SHARES'][r] = shares
        Group_links.rename(columns={'GEOMETRY':'GEOMETRY_ATTRIBUTES'}, inplace = True)
        Street_links = pd.concat([Street_links, Group_links], ignore_index=True)
    else : # non-merged roads
        # Node id generation and save node transform
        [road_node_A,road_node_B] = group_info['NODES_IDS']
        road_node_A = str(road_node_A)
        road_node_B = str(road_node_B)
        node_A_index = nodes_list.index(road_node_A)
        node_B_index = nodes_list.index(road_node_B)
        xa = road_nodes['X'][node_A_index]
        ya = road_nodes['Y'][node_A_index]
        xb = road_nodes['X'][node_B_index]
        yb = road_nodes['Y'][node_B_index]
        s = slope([xa,ya,xb,yb])
        ## add node_A to Groups_nodes
        Group_i_nodes = {}
        Group_i_nodes['GROUP_NODE'] = [nodes_cnt, nodes_cnt + 1]
        Group_i_nodes['X'] = [xa, xb]
        Group_i_nodes['Y'] = [ya, yb]
        Group_i_nodes['GROUP'] = [int(group_i), int(group_i)]
        Group_i_nodes['GROUP_SLOPE'] = [s,s]
        Group_i_nodes['ROAD_NODES'] = [[road_node_A], [road_node_B]]
        Group_i_nodes = pd.DataFrame.from_dict(Group_i_nodes)
        Groups_nodes = pd.concat([Groups_nodes, Group_i_nodes], ignore_index=True)
        ## nodes tranform
        group_nodes_transform[road_node_A]['GROUP_NODES'][group_i] = nodes_cnt
        group_nodes_transform[road_node_B]['GROUP_NODES'][group_i] = nodes_cnt + 1
        group_nodes_transform[road_node_A]['STREET_NODE'] = nodes_cnt
        group_nodes_transform[road_node_B]['STREET_NODE'] = nodes_cnt + 1
        # nodes_transform.loc[node_A_index, 'GROUP_NODE'] = nodes_cnt
        # nodes_transform.loc[node_A_index, 'STREET_NODE'] = nodes_cnt
        # nodes_transform.loc[node_B_index, 'GROUP_NODE'] = nodes_cnt + 1
        # nodes_transform.loc[node_B_index, 'STREET_NODE'] = nodes_cnt + 1
        # Road id generation and save link transform
        ## link in Street_links
        new_link = pd.DataFrame(columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY_ATTRIBUTES'])
        new_link.loc[0, 'STREET_ID'] = links_cnt
        new_link.loc[0, 'NODE_A'] = nodes_cnt
        new_link.loc[0, 'NODE_B'] = nodes_cnt + 1
        new_link.loc[0, 'GEOMETRY_ATTRIBUTES'] = LineString([[xa,ya],[xb,yb]])
        Street_links = pd.concat([Street_links, new_link], ignore_index=True)
        ## link transform
        road_r = group_roads[0]
        r = road_ids.index(road_r)
        links_transform['STREET_IDS'][r] = [links_cnt]
        links_transform['SHARES'][r] = [1]
        links_cnt = links_cnt + 1
        nodes_cnt = nodes_cnt + 2
print(Road_nodes.loc[0,:])

# Step 2) Calculate coordinates of the groups summits at the intersection of streets using Groups_nodes  and replace group_node by street_node in transforms
Intersection_nodes = Road_nodes[(Road_nodes['N_GROUPS'] > 1)]
street_links_node_A = list(Street_links['NODE_A'])
street_links_node_B = list(Street_links['NODE_B'])
Groups_nodes = Groups_nodes.rename(columns={'GROUP_NODE':'NODE_ID'})
# group_node_col = ['GROUP_NODE', 'X', 'Y','GROUP','GROUP_SLOPE','ROAD_NODES']
group_nodes_list = Groups_nodes['NODE_ID']
print(Road_nodes)
print(Intersection_nodes)
print(Groups_nodes)
for road_node_i in Intersection_nodes['NODE_ID'] :
    road_node_i_index = nodes_list.index(road_node_i)
    # 1) Calculate the coordinates of the intersection of the groups/streets
    ## a) get information on group nodes corresponding to the road_node
    Groups_info = pd.DataFrame(columns = Groups_nodes.columns)
    node_i_groups = road_nodes['GROUPS'][road_node_i_index]
    for group_i in node_i_groups :
        subdata = Groups_nodes[Groups_nodes.GROUP == int(group_i)].reset_index(drop=True)
        n_nodes_to_test = subdata.shape[0]
        for j in range(n_nodes_to_test) :
            if road_node_i in list(subdata.loc[j, 'ROAD_NODES']) :
                new_row_i = subdata[subdata.index == j].reset_index(drop=True)
                group_node_i = new_row_i.loc[0,'NODE_ID']
                if group_node_i not in list(Groups_info['NODE_ID']):
                    Groups_info = pd.concat([Groups_info, new_row_i], ignore_index=True)
                other_road_nodes_projected = new_row_i.loc[0, 'ROAD_NODES']
                #other_road_nodes_projected.remove(road_node_i)
                if len(other_road_nodes_projected) > 1 :
                    for road_node_j in other_road_nodes_projected :
                        if road_node_j != road_node_i :
                            road_node_j_index = nodes_list.index(road_node_j)
                            node_j_groups = road_nodes['GROUPS'][road_node_j_index]
                            for group_j in node_j_groups:
                                subdata_j = Groups_nodes[Groups_nodes.GROUP == int(group_j)].reset_index(drop=True)
                                n_nodes_to_test_j = subdata_j.shape[0]
                                for k in range(n_nodes_to_test_j):
                                    if road_node_j in list(subdata_j.loc[k, 'ROAD_NODES']):
                                        new_row_j = subdata_j[subdata_j.index == k].reset_index(drop=True)
                                        group_node_j = new_row_j.loc[0, 'NODE_ID']
                                        if group_node_j not in list(Groups_info['NODE_ID']):
                                            Groups_info = pd.concat([Groups_info, new_row_j], ignore_index=True)
    ## b) get the barycenter of the transformed group nodes
    n_groups = len(Groups_info)
    intersection_x = []
    intersection_y = []
    intersection_coord = []
    for i in range(n_groups) :
        x = Groups_info.loc[i,'X']
        y = Groups_info.loc[i,'Y']
        if (x not in intersection_x) or (y not in intersection_y) :
            intersection_coord.append([x, y])
            intersection_x.append(x)
            intersection_y.append(y)
    points_centered, centroid = Points(intersection_coord).mean_center(return_centroid=True)
    xi = centroid[0]
    yi = centroid[1]
    # 2) Modify the nodes information (id and coordinates) where they appear
    street_node = Groups_info.loc[0,'NODE_ID'] # select the first merging group_node: we will take its id as street id
    for group_node in list(Groups_info['NODE_ID']) :
        ## a) first, in the global list of nodes
        group_node_index = list(Groups_nodes['NODE_ID']).index(group_node)
        Groups_nodes.loc[group_node_index,'NODE_ID'] = street_node # then, the column 'GROUP_NODE' will be labelled 'NODE_ID', refering to the street node id
        Groups_nodes.loc[group_node_index, 'X'] = xi
        Groups_nodes.loc[group_node_index, 'Y'] = yi
        ## b) add id to nodes_transform
        nodes_transform.loc[road_node_i_index, 'STREET_NODE'] = street_node # remark : this is the only case where the street_node ID is different from the group_node_ID
        # group_nodes_transform[road_node_A]['STREET_NODE'] = nodes_cnt
        ## c) modify the node is where it appears in Street_links
        nodes_A_to_modify = get_indexes(street_links_node_A,group_node)
        if len(nodes_A_to_modify) > 0:
            Street_links.loc[nodes_A_to_modify,'NODE_A'] = street_node
        nodes_B_to_modify = get_indexes(street_links_node_B, group_node)
        if len(nodes_B_to_modify) > 0:
            Street_links.loc[nodes_B_to_modify, 'NODE_B'] = street_node

# group_nodes_transform[road_node_A]['GROUP_NODES'] = {g_i : n_i, g_j : n_j}

# Step 3) Generate street links geometry
Street_nodes = Groups_nodes[['NODE_ID', 'X', 'Y']]
n_links = Street_links.shape[0]
for i in range(n_links):
    # for each link, get its nodes ids
    node_A_id = Street_links.loc[i,'NODE_A']
    node_B_id = Street_links.loc[i,'NODE_B']
    # then, their geographical position
    node_A = Street_nodes[Street_nodes.NODE_ID == node_A_id].reset_index()
    node_B = Street_nodes[Street_nodes.NODE_ID == node_B_id].reset_index()
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
print(Street_links['GEOMETRY_ATTRIBUTES'])

# Save files as csv
Street_nodes.to_csv(street_nodes_f, index=False)
Street_links.to_csv(street_links_f, index=False)
pickle.dump(Street_links, open(street_links_dat_f,'wb'))
pickle.dump(links_transform, open(links_transform_f,'wb'))
pickle.dump(nodes_transform, open(nodes_transform_f,'wb'))
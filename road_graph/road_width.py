##################################################################################
# Study: Street network generation
# Purpose: Get mean road width for each road within a cell
# Author: Marjolaine Lannes, adapted from Myrto Valari
# Creation date: January 30, 2023
# Note: data matching apparoach to get the width of OSM roads
##################################################################################
from shapely.geometry import Polygon,Point,LineString
import sys
# sys.path.append("../../../include")
import pandas as pd
import numpy as np
from include.create_poly_around_road import create_1_poly_around_road
from include.find_all_intersections import find_all_intersections
from include.create_circles import create_circles_around_road_nodes
import pickle

# Input and output path
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
roads_bdtopo_dir = path + "temp/load_BDTOPO/subdomains_roads_with_buffer_100m_lamb93/"
roads_OSM_dir = path + "temp/load_OSM/subdomains_roads_lamb93/"
output_dir = path + 'temp/road_graph/subdomains_roads_widths/'

# Get cell information
grid=pickle.load(open(gridfile,'rb'))
poly_num = int(sys.argv[1]) # 4672
x0 = grid['we_id'][poly_num] # index out of range !
y0 = grid['sn_id'][poly_num]

# Input and output files
roads_OSM_f = roads_OSM_dir + 'sub_we_id_from_' + str(x0) + '_sn_id_from_' + str(y0) + '.dat'
roads_bdtopo_f = roads_bdtopo_dir + 'sub_we_id_from0_' + str(x0) + '_sn_id_from0_' + str(y0) + '.dat'
foutput = output_dir + 'sub_we_id_from_' + str(x0) + '_sn_id_from_' + str(y0) + '.dat'

# Load data
roads_OSM_data = pickle.load(open(roads_OSM_f,'rb'))
roads_OSM = pd.DataFrame(roads_OSM_data)
roads_bdtopo_data = pickle.load(open(roads_bdtopo_f,'rb'))
roads_bdtopo = pd.DataFrame(roads_bdtopo_data)
n_roads = roads_OSM.shape[0]

# Median values in BD TOPO roads
median_width = 3.0
median_n_lanes = 2

# Initialize output dict
outdata = {"ROAD_ID": [], "bdtopo_ids": [], "widths": [], "mean_width": []}
for road_id in roads_OSM['ROAD_ID']:
    outdata["ROAD_ID"].append(road_id)
    outdata["bdtopo_ids"].append([])
    outdata["widths"].append([])
    outdata["mean_width"].append(None)

# Data matching algorithm: OSM and BDTOPO roads
for i in range (0,n_roads):
    found = False
    road_xa = roads_OSM.loc[i,'XA']
    road_xb = roads_OSM.loc[i, 'XB']
    road_ya = roads_OSM.loc[i, 'YA']
    road_yb = roads_OSM.loc[i, 'YB']
    road_length = roads_OSM.loc[i, 'LENGTH']
    # Test 1: quasi-identical roads
    radius = 5
    circle_A, circle_B = create_circles_around_road_nodes(road_xa, road_xb, road_ya, road_yb, radius)
    for j, road_bdtopo_id in enumerate(roads_bdtopo['ROAD_ID']) :
        xc = roads_bdtopo.loc[j,'XA']
        yc = roads_bdtopo.loc[j,'YA']
        node_C = Point(xc,yc)
        xd = roads_bdtopo.loc[j,'XB']
        yd = roads_bdtopo.loc[j, 'YB']
        node_D = Point(xd,yd)
        if node_C.intersects(circle_A) and node_D.intersects(circle_B) :
            outdata["bdtopo_ids"][i].append(road_bdtopo_id)
            width = roads_bdtopo.loc[j, 'WIDTH']
            outdata["widths"][i].append(width)
            found = True
        elif node_D.intersects(circle_A) and node_C.intersects(circle_B):
            bd_id = roads_bdtopo.loc[j, 'ROAD_ID']
            outdata["bdtopo_ids"][i].append(bd_id)
            width = roads_bdtopo.loc[j, 'WIDTH']
            outdata["widths"][i].append(width)
            found = True
    ## Then, try a 15m large polygon with 90% coverage
    poly = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, 7.5)  # width = 15 meters
    ids_list, intersections_list, lengths_list,widths_list = find_all_intersections(roads_bdtopo, poly)
    for i_bd, inter in enumerate(intersections_list):
        bd_length = lengths_list[i_bd]
        if inter > 0.9 * bd_length :
            found = True
            outdata["bdtopo_ids"][i].append(ids_list[i_bd])
            outdata["widths"][i].append(widths_list[i_bd])
    # Test 2: nearby roads
    if not found :
        # First try a 25m large polygon with 90% coverage, then 50% coverage
        poly = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, 12.5)  # width = 25 meters
        ids_list, intersections_list, lengths_list,widths_list = find_all_intersections(roads_bdtopo, poly)
        for i_bd, inter in enumerate(intersections_list):
            bd_length = lengths_list[i_bd]
            if inter > 0.9 * bd_length:
                found = True
                outdata["bdtopo_ids"][i].append(ids_list[i_bd])
                outdata["widths"][i].append(widths_list[i_bd])
        if not found :
            for i_bd, inter in enumerate(intersections_list):
                bd_length = lengths_list[i_bd]
                if inter > 0.5 * bd_length:
                    found = True
                    outdata["bdtopo_ids"][i].append(ids_list[i_bd])
                    outdata["widths"][i].append(widths_list[i_bd])
    if not found :
        # Then try a 50m large polygon with 90% coverage, then 50% coverage
        poly = create_1_poly_around_road(road_xa, road_ya, road_xb, road_yb, 25)  # width = 50 meters
        ids_list, intersections_list, lengths_list,widths_list = find_all_intersections(roads_bdtopo, poly)
        for i_bd, inter in enumerate(intersections_list):
            bd_length = lengths_list[i_bd]
            if inter > 0.9 * bd_length:
                found = True
                outdata["bdtopo_ids"][i].append(ids_list[i_bd])
                outdata["widths"][i].append(widths_list[i_bd])
        if not found:
            for i_bd, inter in enumerate(intersections_list):
                bd_length = lengths_list[i_bd]
                if inter > 0.5 * bd_length:
                    found = True
                    outdata["bdtopo_ids"][i].append(ids_list[i_bd])
                    outdata["widths"][i].append(widths_list[i_bd])
    # Case 3: set arbitrary width
    if not found:
        outdata["bdtopo_ids"][i] = None
        outdata["widths"][i].append(median_width*median_n_lanes)
    # Set final width
    outdata["mean_width"][i] = np.mean(outdata["widths"][i])

#Save data
pickle.dump(outdata, open(foutput, 'wb'))
##################################################################################
# Study: Street network generation
# Purpose: Get mean street height for the edges of each street within a cell
# Author: Marjolaine Lannes
# Creation date: January 30
# Note: mean height calculated for each edge of the roads in the cell
# Warning : we do not consider buildings higher than 100 meters
##################################################################################
import sys
import pandas as pd
sys.path.append("../../../include")
from shapely.geometry import Polygon,Point,LineString
from shapely.ops import nearest_points
import numpy as np
from include.create_multiple_poly_around_road import create_left_right_polys_around_road
import pickle

# Load grid data
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
gridfile = path + 'data/grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
grid=pickle.load(open(gridfile,'rb'))

# Get cell information
poly_num = 4672 #int(sys.argv[1])
x0 = grid['we_id'][poly_num]
y0 = grid['sn_id'][poly_num]
prefix = 'we_id_from_' + str(x0) + '_sn_id_from_' + str(y0)

# Input and output files
street_links_f = path + "temp/street_graph/graph_transform/street_links.dat"
streets_f = path + "temp/street_graph/graph_transform/street_links_width_distrib.dat"
buildings_dir = path + 'temp/load_BDTOPO/subdomains_buildings_with_buffer_100m_lamb93/'
fname = 'sub_' + prefix + '.dat'
output_f = path + 'temp/street_graph/attributes/streets_height_width.dat'

# Load data
streets = pickle.load(open(streets_f, 'rb'))
streets_links = pickle.load(open(street_links_f, 'rb'))
buildings = pickle.load(open(buildings_dir + fname, 'rb'))
streets_list = list(streets['STREET_ID'])
n_streets = len(streets_list)
num_build=len(buildings['BUILDING_ID'])
print('This processor will loop over '+ str(num_build) + ' buildings')
print('This processor will loop over ' + str(n_streets) + ' streets.')

# Calculate street heights in BD TOPO
print(streets.columns)
init_list = [0.0]*n_streets
streets_height_width = pd.DataFrame(columns=['STREET_ID', 'HEIGHT_LEFT', 'HEIGHT_RIGHT', 'WIDTH'])
streets_height_width['STREET_ID'] = streets_list
for col in ['HEIGHT_LEFT', 'HEIGHT_RIGHT', 'WIDTH']:
    streets_height_width[col] = init_list
# streets_height_width = {'STREET_ID':streets_list, 'HEIGHT_LEFT': init_list, 'HEIGHT_RIGHT':init_list, 'WIDTH':init_list}
for s, street_id in enumerate(streets_list) :
    print("street ",s)
    width = streets.loc[s,'TOTAL_WIDTH']
    street_line = streets_links.loc[s,'GEOMETRY_ATTRIBUTES']
    street_line_str = str(street_line)
    [xa, ya], [xb,yb] = street_line.coords
    # xa = streets.loc[s, 'XA']
    # xb = streets.loc[s, 'XB']
    # ya = streets.loc[s, 'YA']
    # yb = streets.loc[s, 'YB']
    # create a polygon around road, including bicycle paths and surrounding buildings
    dist = width / 2 + 10  # (meters) 10 meters is the building width
    [road_polys_left, road_polys_right, n_poly, road_length] = create_left_right_polys_around_road(xa, ya, xb, yb, dist)
    heights_left = [0.0] * n_poly
    n_buildings_left = [0] * n_poly
    width_left = np.nan
    heights_right = [0.0] * n_poly
    n_buildings_right = [0] * n_poly
    width_right = np.nan
    # get buildings height inside the road polygon for each side of the road
    # print(buildings['POLYGON'])
    for bld_cnt, building_poly in enumerate(buildings['POLYGON']):
        for i, road_poly in enumerate(road_polys_right):
            if building_poly.intersects(road_poly):
                h = buildings['HEIGHT'][bld_cnt]
                if not (np.isnan(h) or h > 100):
                    heights_right[i] += float(h)
                    n_buildings_right[i] += 1
                    [P1, P2] = list(nearest_points(building_poly, street_line))
                    new_width = P1.distance(P2)
                    width_right = np.nanmin([width_right, new_width])
        for i, road_poly in enumerate(road_polys_left):
            if building_poly.intersects(road_poly):
                h = buildings['HEIGHT'][bld_cnt]
                if not np.isnan(h) or h > 100:
                    heights_left[i] += float(h)
                    n_buildings_left[i] += 1
                    [P1, P2] = list(nearest_points(building_poly, street_line))
                    new_width = P1.distance(P2)
                    width_left = np.nanmin([width_left, new_width])
    # average the heights if there are some buildings
    n_buildings = sum(n_buildings_left) + sum(n_buildings_right)
    if n_buildings > 0:
        # average height on the left side of the road
        aver_list_left = [0.0] * n_poly
        for j in range(0, n_poly):
            if n_buildings_left[j] > 0:
                aver_list_left[j] = heights_left[j] / n_buildings_left[j]
            else:
                aver_list_left[j] = 0
        aver_left = sum(aver_list_left) / n_poly
        # average height on the right side of the road
        aver_list_right = [0.0] * n_poly
        for k in range(0, n_poly):
            if n_buildings_right[k] > 0:
                aver_list_right[k] = heights_right[k] / n_buildings_right[k]
            else:
                aver_list_right[k] = 0
        aver_right = sum(aver_list_right) / n_poly
    # If no building found, try with a larger polygon around the road
    else:
        # print('No bulding found for street id '+str(road_id))
        dist = width / 2 + 20  # (meters) 20 meters is the building width
        [road_polys_left, road_polys_right, n_poly, road_length] = create_left_right_polys_around_road(xa, ya, xb, yb,
                                                                                                       dist)
        heights_left = [0.0] * n_poly
        n_buildings_left = [0] * n_poly
        heights_right = [0.0] * n_poly
        n_buildings_right = [0] * n_poly
        # get buildings height inside the road polygon for each side of the road
        for bld_cnt, building_poly in enumerate(buildings['POLYGON']):
            for i, road_poly in enumerate(road_polys_right):
                if building_poly.intersects(road_poly):
                    h = buildings['HEIGHT'][bld_cnt]
                    if not np.isnan(h) or h > 100:
                        heights_right[i] += float(h)
                        n_buildings_right[i] += 1
                        [P1, P2] = list(nearest_points(building_poly, street_line))
                        new_width = P1.distance(P2)
                        width_right = np.nanmin([width_right, new_width])
            for i, road_poly in enumerate(road_polys_left):
                if building_poly.intersects(road_poly):
                    h = buildings['HEIGHT'][bld_cnt]
                    if not np.isnan(h) or h > 100:
                        heights_left[i] += float(h)
                        n_buildings_left[i] += 1
                        [P1, P2] = list(nearest_points(building_poly, street_line))
                        new_width = P1.distance(P2)
                        width_left = np.nanmin([width_left, new_width])
        # average if there are some buildings
        n_buildings = sum(n_buildings_left) + sum(n_buildings_right)
        if n_buildings > 0:
            # average height on the left side of the road
            aver_list_left = [0.0] * n_poly
            for i in range(0, n_poly):
                if n_buildings_left[i] > 0:
                    aver_list_left[i] = heights_left[i] / n_buildings_left[i]
                else:
                    aver_list_left[i] = 0
            aver_left = sum(aver_list_left) / n_poly
            # average height on the right side of the road
            aver_list_right = [0.0] * n_poly
            for i in range(0, n_poly):
                if n_buildings_right[i] > 0:
                    aver_list_right[i] = heights_right[i] / n_buildings_right[i]
                else:
                    aver_list_right[i] = 0
            aver_right = sum(aver_list_right) / n_poly
        # If no building found (again), then set height = 0
        else:
            aver_left = 0.0
            aver_right = 0.0
    # Save the collected data in the output dataframe
    streets_height_width.loc[s,'HEIGHT_LEFT'] = aver_left
    streets_height_width.loc[s,'HEIGHT_RIGHT'] = aver_right
    if np.isnan(width_left) or np.isnan(width_right) :
        streets_height_width.loc[s,'WIDTH'] = np.nanmax([width_left, width_right])
    else :
        streets_height_width.loc[s,'WIDTH'] = width_left + width_right
print(streets_height_width)

# Save height and width data for each road from BD TOPO (IDF)
fout = open(output_f, 'wb')
pickle.dump(streets_height_width, fout)
fout.close()

##################################################################################
# Study: Street network generation
# Purpose: Get distance to buildings and mean street height for the edges of each street within a cell
# Note: mean height and distance to buildings are calculated for left and right edges of the roads in the cell
# Warning : we do not consider buildings higher than 100 meters
##################################################################################
import sys
import pandas as pd
import geopandas as gpd
sys.path.append("../../include")
from shapely.geometry import Polygon,Point,LineString
from shapely.ops import nearest_points
from shapely.wkt import loads
import numpy as np
from create_multiple_poly_around_road import create_left_right_polys_around_road
from create_poly_around_road import create_1_poly_around_road
from get_indexes import get_indexes
import pickle, os, time
from ast import literal_eval
start_time = time.time()

# Parameters
max_height = 100 #meters
max_width_then_open = 160 #meters
street_crossed_by_buidings_ratio = 0.5 # if more than n% of the street geometry crosses the building, take the main road geometry instead

# Directories
path = "../../"
input_dir = path + "data/"
cache_dir  = path + "temp/"

# Get cell information
gridfile = input_dir + 'grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
grid=pickle.load(open(gridfile,'rb'))
poly_num = int(sys.argv[1])
x0 = grid['we_id'][poly_num]
y0 = grid['sn_id'][poly_num]

# Input and output files
buildings_dir = cache_dir + 'data/urban_form/subdomains_buildings_with_buffer_100m_lamb93/'
groups_f = cache_dir + "street_graph/groups/groups.dat"
links_inverse_transform_f = cache_dir + "street_graph/graph_transform/links_inverse_transform.csv"
road_links_f = cache_dir + "data/load_OSM/roads_links.csv"
streets_dir = cache_dir + "street_graph/graph_transform/street_cells/"
fname = 'sub_we_id_from_' + str(x0) + '_sn_id_from_' + str(y0) + '.dat'
fnameb = 'sub_we_id_from0_' + str(x0) + '_sn_id_from0_' + str(y0) +'.dat'
outdir = cache_dir + 'street_graph/attributes/subdomains_streets_height_width/'

# Load data
streets = pickle.load(open(streets_dir + fname, 'rb'))
streets_list = list(streets['STREET_ID'])
n_streets = len(streets_list)
if n_streets > 0 :
    buildings = pickle.load(open(buildings_dir + fnameb, 'rb'))
    buildings_df = pd.DataFrame(buildings)
    buildings_gdf = gpd.GeoDataFrame(buildings_df, crs='epsg:2154', geometry='POLYGON')

# Calculate street heights in BD TOPO
init_list = [0]*n_streets
streets_height_width = {'STREET_ID':streets_list, 'HEIGHT_LEFT': init_list.copy(), 'HEIGHT_RIGHT':init_list.copy(), 'WIDTH':init_list.copy(),
                        'OPEN':['no']*n_streets, 'RATIO': init_list.copy()}
width_threshold = max_width_then_open / 2
for s, street_id in enumerate(streets_list) :
    print(s)
    is_tunnel = streets['TUNNEL'][s]
    if is_tunnel :
        streets_height_width['HEIGHT_LEFT'][s] = np.nan
        streets_height_width['HEIGHT_RIGHT'][s] = np.nan
        streets_height_width['WIDTH'][s] = np.nan
    else :
        xa = streets['XA'][s]
        ya = streets['YA'][s]
        xb = streets['XB'][s]
        yb = streets['YB'][s]
        street_line = streets['GEOMETRY'][s]
        coverage_left = [0] # i-th value = 1 if a building has been found in the i-th left polygon, 0 otherwise
        coverage_right = [0] # idem right
        n_poly = 10
        mean_width = np.nan
        building_distance = 10 #(meters) for the building width
        while (building_distance < width_threshold) & (sum(coverage_left) < 0.8 * n_poly) & (sum(coverage_right) < 0.8 * n_poly) :
            # create a polygon around road, including bicycle paths and surrounding buildings
            [road_polys_left, road_polys_right, centers_update, n_poly, road_length] = create_left_right_polys_around_road(xa,ya,xb,yb,building_distance)
            street_poly = create_1_poly_around_road(xa, ya, xb, yb, building_distance)
            SUBbuildings = buildings_gdf[buildings_gdf.POLYGON.intersects(street_poly)]
            if building_distance == 10 : # mark for which polygons the width/height were found
                coverage_right = [0] * n_poly
                coverage_left = [0] * n_poly
                heights_left = [0.0] * n_poly
                n_buildings_left = [0] * n_poly
                widths_left = [np.nan] * n_poly
                heights_right = [0.0] * n_poly
                n_buildings_right = [0] * n_poly
                widths_right = [np.nan] * n_poly
                centers = centers_update
                for i in range(n_poly):
                    if SUBbuildings[SUBbuildings.POLYGON.intersects(centers[i])].shape[0] > 0:
                        coverage_right[i] = 1
                        coverage_left[i] = 1
                # if more than 50% of the street geometry crosses the building, take the main road geometry instead
                if sum(coverage_right) / n_poly > street_crossed_by_buidings_ratio :
                    if not "roads" in globals():
                        roads = pd.read_csv(road_links_f).astype({'ID': str})
                        roads['GEOMETRY'] = roads['GEOMETRY'].apply(lambda x: loads(x))
                        roads_list = list(roads['ID'])
                        Groups = pickle.load(open(groups_f, 'rb'))
                        links_inverse_transform = pd.read_csv(links_inverse_transform_f)
                        links_inverse_transform['GROUP_ID'] = links_inverse_transform['GROUP_ID'].apply(
                            lambda x: literal_eval(x))
                    links_inverse_group_g = links_inverse_transform[links_inverse_transform.STREET_ID == street_id].reset_index(drop=True)
                    if len(links_inverse_group_g.loc[0, 'GROUP_ID']) > 0 :
                        group_g = links_inverse_group_g.loc[0, 'GROUP_ID'][0]
                        matching_road = Groups[group_g]['MAIN_ROAD_IDS'][0]
                        road_index = roads_list.index(str(matching_road))
                        xa = roads.loc[road_index, 'XA']
                        ya = roads.loc[road_index, 'YA']
                        xb = roads.loc[road_index, 'XB']
                        yb = roads.loc[road_index, 'YB']
                        line = roads.loc[road_index, 'GEOMETRY']
                        n_poly_old = n_poly
                        [road_polys_left, road_polys_right, centers_update, n_poly, road_length] = create_left_right_polys_around_road(xa, ya, xb, yb, building_distance)
                        points = [line.interpolate((i / n_poly), normalized=True) for i in range(0, n_poly + 1)]
                        centers = [LineString([points[i], points[i+1]]).centroid for i in range(0, n_poly)]
                        street_poly = create_1_poly_around_road(xa, ya, xb, yb, building_distance)
                        SUBbuildings = buildings_gdf[buildings_gdf.POLYGON.intersects(street_poly)]
                        if n_poly != n_poly_old :
                            coverage_right = [0] * n_poly
                            coverage_left = [0] * n_poly
                            heights_left = [0.0] * n_poly
                            n_buildings_left = [0] * n_poly
                            widths_left = [np.nan] * n_poly
                            heights_right = [0.0] * n_poly
                            n_buildings_right = [0] * n_poly
                            widths_right = [np.nan] * n_poly
                            centers = centers_update
                            for i in range(n_poly):
                                if SUBbuildings[SUBbuildings.POLYGON.intersects(centers[i])].shape[0] > 0:
                                    coverage_right[i] = 1
                                    coverage_left[i] = 1
            # get buildings height/width inside the road polygon for each side of the road
            if SUBbuildings.shape[0] > 0:
                for i in get_indexes(coverage_right,0):
                    road_poly = road_polys_right[i]
                    intersected_buildings = SUBbuildings[SUBbuildings.POLYGON.intersects(road_poly)].reset_index(drop=True)
                    new_width = widths_right[i]
                    for bld_cnt, building_poly in enumerate(intersected_buildings['POLYGON']):
                        h = intersected_buildings.loc[bld_cnt,'HEIGHT']
                        if not (np.isnan(h) or h > max_height):
                            heights_right[i] += float(h)
                            n_buildings_right[i] += 1
                            coverage_right[i] = 1
                            [P1, P2] = list(nearest_points(building_poly, centers[i]))
                            new_width = np.nanmin([new_width, P1.distance(P2)])
                    widths_right[i] = new_width
                for i in get_indexes(coverage_left,0):
                    road_poly = road_polys_left[i]
                    intersected_buildings = SUBbuildings[SUBbuildings.POLYGON.intersects(road_poly)].reset_index(drop=True)
                    new_width = widths_left[i]
                    for bld_cnt, building_poly in enumerate(intersected_buildings['POLYGON']):
                        h = intersected_buildings.loc[bld_cnt, 'HEIGHT']
                        if not (np.isnan(h) or h > max_height):
                            heights_left[i] += float(h)
                            n_buildings_left[i] += 1
                            coverage_left[i] = 1
                            [P1, P2] = list(nearest_points(building_poly, centers[i]))
                            new_width = np.nanmin([new_width, P1.distance(P2)])
                    widths_left[i] = new_width
            building_distance += 10 # If no building found, try with a larger polygon around the road
        if sum(n_buildings_left) + sum(n_buildings_right) > 0:
            # average height on the left side of the road
            aver_list_left = [0.0] * n_poly
            for j in range(0, n_poly):
                if n_buildings_left[j] > 0: # IndexError: list index out of range
                    aver_list_left[j] = heights_left[j] / n_buildings_left[j]
                else:
                    aver_list_left[j] = 0
            aver_height_left = sum(aver_list_left) / n_poly
            # average height on the right side of the road
            aver_list_right = [0.0] * n_poly
            for k in range(0, n_poly):
                if n_buildings_right[k] > 0:
                    aver_list_right[k] = heights_right[k] / n_buildings_right[k]
                else:
                    aver_list_right[k] = 0
            aver_height_right = sum(aver_list_right) / n_poly
            # calculate road width
            if (aver_height_right > 0) and (aver_height_left > 0) :
                widths = [sum(x) for x in zip(widths_left, widths_right)]
                mean_width = np.nanmedian(widths)
            elif aver_height_right == 0.0:
                mean_width = 2 * np.nanmedian(widths_left)
            elif aver_height_left == 0.0:
                mean_width = 2 * np.nanmedian(widths_right)
        else: # If no building found (again), then set height = 0
            aver_height_left = 0.0
            aver_height_right = 0.0
            mean_width = np.nan
            # Save the collected data in the output dataframe
        streets_height_width['HEIGHT_LEFT'][s] = aver_height_left
        streets_height_width['HEIGHT_RIGHT'][s] = aver_height_right
        streets_height_width['WIDTH'][s] = mean_width

# Save height and width data for each road
fout = open(outdir + fname, 'wb')
pickle.dump(streets_height_width, fout)
fout.close()

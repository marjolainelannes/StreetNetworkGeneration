##################################################################################
# Study: Street network generation
# Purpose: For each road link, save list of links that are close enough
# Notes: (*) close in the sense of to a certain number of nodes (not metric)
# (**) They are used to accelerate road groups identification easy to parallelize if necessary
##################################################################################
import pandas as pd
import geopandas as gpd
from shapely.wkt import loads
from shapely.geometry import *
import pickle
import time

# Config
path = "../../../"
input_dir = path + "data/"
cache_dir  = path + "temp/"
links_distance = 500 # meters

# Input/output files
links_f = cache_dir + "data/load_OSM/roads_links.csv"
output_connections_f = cache_dir + "data/load_OSM/neighborhood_links_" + str(links_distance)+ "m.dat"

# Load network information from OSM
Links = pd.read_csv(links_f) # link_col = ['ID', 'NODE_A', 'NODE_B', 'XA', 'YA', 'XB', 'YB', 'LENGTH', 'ROUNDABOUT', 'TUNNEL', 'LINE']
Links = Links[['ID', 'LINE']]
Links = Links.astype({'LINE':'str', 'ID':'str'})
Links.LINE = Links.LINE.apply(loads)
Links_gdf = gpd.GeoDataFrame(Links,crs='epsg:2154', geometry = 'LINE')
print(Links_gdf)
print("number of links:", Links_gdf.shape[0])

# Get list of close links
start_time = time.time()
connected_links = {}
for i, link_geom in enumerate(Links_gdf['LINE']) :
    link_id = Links_gdf.loc[i,'ID']
    center = link_geom.centroid
    buffer_radius = link_geom.length + links_distance
    poly = center.buffer(buffer_radius)
    subdata = Links_gdf[Links_gdf.geometry.intersects(poly)]
    connected_links[link_id] = list(subdata.ID)
    if i % 10000 == 0 :
        print(i)
        time_one_cell = time.gmtime(time.time() - start_time)
        print("It took", time_one_cell.tm_hour, "h", time_one_cell.tm_min, "min", time_one_cell.tm_sec, "seconds")

# Save data
pickle.dump(connected_links, open(output_connections_f, 'wb'))

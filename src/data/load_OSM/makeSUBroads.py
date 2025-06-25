##################################################################################
# Study: Street network generation
# Purpose: subdivide roads information by cells for parallelization of road width calculus
# Note: Get roads information (id, geometry) of OSM roads for each cell of the grid
##################################################################################
import pandas as pd
import pickle, time
import geopandas as gpd
from shapely.wkt import loads

# Param
path = "../../../"
temp_dir = path + "temp/"
data_dir = path + "data/"

# Input/output files
roads_f = temp_dir + "data/load_OSM/roads_links.csv"
grid_file = data_dir + 'grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
cells_roads_f = path + "temp_VF/data/load_OSM/roads_list_per_cell.dat"

# Load data
roads_df = pd.read_csv(roads_f)
roads_df = roads_df.astype({'GEOMETRY':'str', 'ID':'str'})
roads_df.GEOMETRY = roads_df.GEOMETRY.apply(loads)
roads = gpd.GeoDataFrame(roads_df,crs='epsg:2154', geometry = 'GEOMETRY')
n_roads = roads.shape[0]
print("Number of roads: ",n_roads)

# Save data for the roads within each cell
start = time.time()
grid = pickle.load(open(grid_file, 'rb'))
cells_roads_list = {}
for poly_num, poly in enumerate(grid['polygon']):
    # Select roads intersecting the cell
    subdata = roads[roads.geometry.intersects(poly)]
    cells_roads_list[poly_num] = list(subdata.ID)
    if poly_num % 100 == 0:
        timer = time.gmtime(time.time() - start)
        print("time for ", poly_num, " cells:")
        print("It took", timer.tm_hour, "h", timer.tm_min, "min", timer.tm_sec, "seconds")

# Save the lists of cells for each road
pickle.dump(cells_roads_list, open(cells_roads_f, 'wb'))
print('Cells added to road data...')

##################################################################################
# Study: Exposure modelling framework
# Purpose: to save buildings information from BDNB
# Note: Get buildings information (id, height, geometry) for each cell of the grid
##################################################################################

from shapely.geometry import Polygon,Point,LineString
from shapely.wkt import loads
import pickle, os
import pandas as pd
import geopandas as gpd
import time
start_time = time.time()

# Config
path = "../../../"
input_dir = path + "data/"
cache_dir  = path + "temp/"

# Input and output files
buildings_input = open(cache_dir + 'data/urban_form/buildings_BDNB.dat','rb')
grid_file = open(input_dir + 'grid_data/polygons_cells_IDF_with_buffer_100m_lambert93.dat','rb')
outdir = cache_dir + 'data/urban_form/subdomains_buildings_with_buffer_100m_lamb93/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Load buildings input
grid = pickle.load(grid_file)
buildings = pickle.load(buildings_input)
buildings_df = pd.DataFrame(buildings)
buildings_df['GEOMETRY'] = buildings_df['GEOMETRY'].apply(loads)
data = gpd.GeoDataFrame(buildings_df, geometry='GEOMETRY')
fields = data.columns
buildings_input.close()
grid.close()
print('Buildings data loaded...')

# Buildings Polygons with 100m buffer
print("Building Polygons with 100m buffer")
for poly_num, poly in enumerate(grid['polygon']) :
    # get cell coordinates
    x = grid['we_id'][poly_num]
    y = grid['sn_id'][poly_num]
    # one file per cell with all buildings included in it
    prefix = 'we_id_from0_' + str(x) + '_sn_id_from0_' + str(y)
    foutname = open(outdir + 'sub_' + prefix + '.dat', 'wb')
    # select all buildings within the cell
    subdata = data[data.geometry.within(poly)]
    print('Buildings data croped in cells...')
    # get their id, height and geometry
    buildings_in_cell = {'BUILDING_ID':list(subdata.ID) , 'HEIGHT':list(subdata.HEIGHT) , 'POLYGON':list(subdata.GEOMETRY) }
    print('Now saving....')
    # save it
    pickle.dump(buildings_in_cell, foutname)
    foutname.close()
    print('Done with cell ', poly_num)

# Time
print('All buildings from BDNB are loaded.')
time_tot_sec = time.time() - start_time
time_tot = time.gmtime(time_tot_sec)
print("Total time: ", time_tot.tm_hour, 'hour', time_tot.tm_min, 'min', time_tot.tm_sec, 'sec.')

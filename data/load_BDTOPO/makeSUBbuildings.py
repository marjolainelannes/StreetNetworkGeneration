##################################################################################
# Study: Exposure modelling framework
# Purpose: to save buildings information from BDTOPO
# Author: Marjolaine Lannes, adapted from Myrto Valari
# Creation date: November 24, 2022
# Note: Get buildings information (id, height, geometry) for each cell of the grid
##################################################################################

from shapely.geometry import Polygon,Point,LineString
import pickle
import geopandas as gpd

# Input and output files
path = "../../../"
buildings_input = path + 'data/BDTOPO/BATIMENT.shp'

# Load buildings input
data = gpd.read_file(buildings_input)
fields = data.columns
print('Buildings data loaded...')

for buffer in ['with_buffer_100m', 'no_buffer']:
    print(buffer)
    # Load grid file
    grid_file = path + 'data/grid_data/polygons_cells_IDF_' + buffer + '_lambert93.dat'
    outdir = path + 'temp/load_BDTOPO/subdomains_buildings_' + buffer + '_lamb93/'
    grid = pickle.load(open(grid_file, 'rb'))
    # Building Polygons
    for poly_num, poly in enumerate(grid['polygon']) :
        # get cell coordinates
        x = grid['we_id'][poly_num]
        y = grid['sn_id'][poly_num]
        # one file per cell with all buildings included in it
        prefix = 'we_id_from0_' + str(x) + '_sn_id_from0_' + str(y)
        foutname = outdir + 'sub_' + prefix + '.dat'
        # select all buildings within the cell
        subdata = data[data.geometry.within(poly)]
        print('Buildings data croped in cells...')
        # get their id, height and geometry
        buildings_in_cell = {'BUILDING_ID':list(subdata.index) , 'HEIGHT':list(subdata.HAUTEUR) , 'POLYGON':list(subdata.geometry) }
        print('Now saving....')
        # save it
        pickle.dump(buildings_in_cell, open(foutname, 'wb'))
        print('Done with cell ', poly_num)

print('All buildings from BD TOPO are loaded.')
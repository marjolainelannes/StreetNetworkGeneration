##################################################################################
# Study: Street network generation
# Purpose: to save buildings information from BDTOPO
# Author: Marjolaine Lannes
# Creation date: June 15, 2023
# Note: Get buildings information (id, height, geometry) for Ile-de-France region
##################################################################################

from shapely.geometry import Polygon,Point,LineString
import pickle
import geopandas as gpd

# Input and output files
path = "../../../"
buildings_input = path + 'data/BDTOPO/BATIMENT.shp'
foutname = path + 'temp/load_BDTOPO/buildings.dat'

# Load buildings input
data = gpd.read_file(buildings_input)
fields = data.columns
print('Buildings data loaded...')

# Save information
buildings = {'BUILDING_ID':list(data.index) , 'HEIGHT':list(data.HAUTEUR) , 'POLYGON':list(data.geometry) }
pickle.dump(buildings, open(foutname, 'wb'))
print('All buildings from BD TOPO are loaded.')
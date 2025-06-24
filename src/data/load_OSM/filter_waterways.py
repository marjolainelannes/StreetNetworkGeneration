##################################################################################
# Study: Street network generation
# Purpose: to save waterways information from OSM
# Author: Marjolaine Lannes
# Creation date: September 12, 2023
# Note: Get (filtered) waterways location for Ile-de-France region
##################################################################################
import pickle
import geopandas as gpd

# Input and output files
path = "../../../"
waterways_input = path + 'data/OSM/waterways/waterways.shp'
fout_dat = path + 'temp/load_OSM/waterways/waterways.dat'
fout_shp = path + 'temp/load_OSM/waterways/filtered_waterways.shp'

# Load buildings input and filter underground / overground / tunnels
data = gpd.read_file(waterways_input)
filtered_waterways = data[~ data['location'].isin(["underground","overground"])]
filtered_waterways = filtered_waterways[~filtered_waterways['man_made'].str.contains("pipeline",na=False)]
filtered_waterways = filtered_waterways[~filtered_waterways['tunnel'].isin(["yes", "culvert"])]
print('Waterways data loaded...')

# Save information
waterways = {'WATERWAY_ID':list(filtered_waterways.full_id), 'GEOMETRY':list(filtered_waterways.geometry) }
pickle.dump(waterways, open(fout_dat, 'wb'))
filtered_waterways.to_file(fout_shp)
print('All waterways from OSM are loaded.')

##################################################################################
# Study: Street network generation
# Author: Marjolaine Lannes
# Creation date: December 12, 2022
##################################################################################
import geopandas as gpd
import pickle,sys

# Input and output files
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
sidewalks_file = path + 'temp/load_sidewalks/sidewalks_lamb93.shp'
foutname = path + 'temp/load_sidewalks/Paris_sidewalks.dat'

# Load sidewalks sidewalks_data
sidewalks_data=gpd.read_file(sidewalks_file)
print('Shapely file loaded...')

outdata={'ID':[],'WIDTH':[],'LINE':[]}
for i,feat in sidewalks_data.iterrows() :
    line=feat.geometry
    outdata['ID'].append(sidewalks_data.ID[i])
    outdata['WIDTH'].append(sidewalks_data.WIDTH[i])
    outdata['LINE'].append(line)
print('Now saving...')
pickle.dump(outdata,open(foutname,'wb'))
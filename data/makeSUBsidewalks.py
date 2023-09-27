##################################################################################
# Study: Exposure modelling framework
# Purpose: For a cell, read the sidewalks and select those in the the 1km cells with 1000m buffer (3km x 3km)
# Author: Myrto Valari
# Creation date: December 12, 2022
##################################################################################

import geopandas as gpd
import pickle,sys

## Only with buffer
buffer = "with_buffer_100m"

# Input and output files
path = "C:/Users/marjolaine.lannes/PycharmProjects/Road_network_topology/"
gridfile = path + 'data/grid_data/polygons_cells_IDF_' + buffer + '_lambert93.dat'
sidewalks_file = path + 'temp/load_sidewalks/sidewalks_lamb93.shp'
outdir = path + 'temp/load_sidewalks/sub_sidewalks_' + buffer + '_lamb93/'

# Load grid data
grid=pickle.load(open(gridfile,'rb'))

# Get the cell information
poly_num = int(sys.argv[1])
poly=grid['polygon'][poly_num]
x=grid['we_id'][poly_num]
y=grid['sn_id'][poly_num]
prefix='we_id_from0_' + str(x) +'_sn_id_from0_' + str(y)
print(poly_num)

# Load sidewalks sidewalks_data
sidewalks_data=gpd.read_file(sidewalks_file)
print('Shapely file loaded on ' + str(poly_num))

# Get & save all sidewalks characteristics for those intersecting the cell
foutname=outdir+'sub_'+prefix+'.dat'
outdata={'ID':[],'HERID':[],'HERID1':[],'WIDTH':[],'LINE':[]}
for i,feat in sidewalks_data.iterrows() :
        line=feat.geometry
        if poly.intersects(line) :
           outdata['ID'].append(sidewalks_data.ID[i])
           outdata['HERID'].append(sidewalks_data.HERID[i])
           outdata['HERID1'].append(sidewalks_data.HERID1[i])
           outdata['WIDTH'].append(sidewalks_data.WIDTH[i])
           outdata['LINE'].append(line)
print('Now saving...' + str(poly_num))
pickle.dump(outdata,open(foutname,'wb'))
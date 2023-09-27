##################################################################################
# Study: Exposure modelling framework
# Purpose: to save roads information from BDTOPO
# Author: Marjolaine Lannes, adapted from Myrto Valari
# Creation date: November 24, 2022
# Note: Get roads information (id, width, lanes, geometry) for one cell of the grid
##################################################################################
from shapely.geometry import Polygon, Point, LineString
import pickle, sys
import numpy as np
import geopandas as gpd
import time

# Input files
path = "C:/Users/marjolaine.lannes/PycharmProjects/Road_network_topology/"
roads_input_file = path + 'data/BDTOPO/TRONCON_DE_ROUTE.shp'
poly_num = int(sys.argv[1])

# Load roads_data
roads_data = gpd.read_file(roads_input_file)
print(roads_input_file + ' loaded')
n_roads = roads_data.shape[0]

# Define road types to collect / exclude
outdata = {'ROAD_ID': [],'IDBD':[], 'WIDTH': [], 'N_LANES': [], 'XA': [], 'YA': [], 'XB': [], 'YB': [], 'LINE': []}
to_exclude = ['Sentier', 'Escalier', 'Route empierrée', 'Piste cyclable', 'Chemin', 'Bac ou liaison maritime']
# remaining road types: ['Route à 1 chaussée', 'Route à 2 chaussées', 'Rond-point', 'Bretelle', 'Type autoroutier']

# Save data for the roads within each cell
start = time.time()
buffer = 'with_buffer_100m'
# Load grid file
grid_file = path + 'data/grid_data/polygons_cells_IDF_' + buffer + '_lambert93.dat'
outdir = 'temp/load_BDTOPO/subdomains_roads_' + buffer + '_lamb93/'
grid = pickle.load(open(grid_file, 'rb'))

# Get the cell information
x = grid['we_id'][poly_num]
y = grid['sn_id'][poly_num]
poly = grid['polygon'][poly_num]
## one file per cell with all roads included in it
prefix = 'we_id_from0_' + str(x) + '_sn_id_from0_' + str(y)
foutname = outdir + 'sub_' + prefix + '.dat'

# Save road network data
myid = 0
for i, road in roads_data.iterrows(): # for each road in BD TOPO database, check if it actually is a road & in the cell
    road_type = road.NATURE
    xis, yis = road.geometry.xy
    if road_type not in to_exclude:
        line = road.geometry
        if poly.intersects(line): # if it is, get its information
            n_lanes = road.NB_VOIES
            width = road.LARGEUR
            idbd = road.ID
            if width is None or np.isnan(width):
                width = -999.0
            if n_lanes is None or n_lanes == 0 :
                n_lanes = -999.0
            if len(xis) == 2:
                x1, x2 = xis[:]
                y1, y2 = yis[:]
                outdata['ROAD_ID'].append(myid)
                outdata['IDBD'].append(idbd)
                outdata['N_LANES'].append(n_lanes)
                outdata['WIDTH'].append(width)
                outdata['XA'].append(x1)
                outdata['YA'].append(y1)
                outdata['XB'].append(x2)
                outdata['YB'].append(y2)
                outdata['LINE'].append(line)
                myid += 1
            else:
                a = list(zip(xis, yis))
                x1, y1 = a[0]
                for j in range(len(xis) - 1):
                    x2, y2 = a[j + 1] ########################## MODIF
                    outdata['ROAD_ID'].append(myid)
                    outdata['IDBD'].append(idbd)
                    outdata['N_LANES'].append(n_lanes)
                    outdata['WIDTH'].append(width)
                    outdata['XA'].append(x1)
                    outdata['YA'].append(y1)
                    outdata['XB'].append(x2)
                    outdata['YB'].append(y2)
                    line = LineString(((x1, y1), (x2, y2)))
                    outdata['LINE'].append(line)
                    myid += 1
                    x1, y1 = x2, y2
        else : # to take into account all the roads parts we have ignored bc outside the cell
            for j in range(len(xis)-1):
                myid += 1
    else : # to take into account all the roads parts we have ignored bc of the road type
        for j in range(len(xis)-1) :
            myid += 1
# Save it
pickle.dump(outdata, open(foutname, 'wb'))
end = time.time()
print('Done with cell ', poly_num, "in ", (end - start), " seconds" )
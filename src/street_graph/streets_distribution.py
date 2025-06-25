##################################################################################
# Study: Street network generation
# Purpose: Get streets attributes: width, height, cells distribution
##################################################################################
import pickle, os
import pandas as pd
import geopandas as gpd
from itertools import repeat
import time
import shapely.ops as sp_ops
import pyproj
from shapely import wkt

# Directories
path = "../../"
data_dir  = path + "data/"
temp_dir  = path + "temp/"

# Input and output files
grid_file = data_dir + 'grid_data/polygons_cells_IDF_no_buffer_lambert93.dat'
transform_dir = temp_dir + "streets_graph/graph_transform/"
street_nodes_f = transform_dir + "street_nodes.csv"
street_links_f = transform_dir + "street_links.csv"
links_transform_f = transform_dir + "links_transform.dat"
links_inverse_transform_f = transform_dir + "links_inverse_transform.dat"
street_links_attributes_f = transform_dir + "street_links_distrib.csv"
street_cells_distrib_f = transform_dir + "street_cells_distrib.csv"
outdir = transform_dir + 'street_cells/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Load data
Street_links = pd.read_csv(street_links_f) # columns=['STREET_ID', 'NODE_A', 'NODE_B', 'GEOMETRY','AXIAL_REPRESENTATION, 'XA', 'XB', 'YA', 'YB']
Street_nodes = pd.read_csv(street_nodes_f) # columns = ['NODE_ID', 'X', 'Y']
links_pkl = open(links_transform_f,'r')
links_transform = pickle.load(open(links_transform_f,'r'))
links_pkl.close() # columns=['ROAD_ID', 'GROUP_ID', 'PARALLEL', 'MERGED', 'COLOR', 'GEOMETRY','STREET_IDS', 'SHARES']
grid = pickle.load(open(grid_file, 'rb'))
my_transformer = pyproj.Transformer.from_crs('EPSG:2154', 'EPSG:3857', always_xy=True)

# 1) Inverse transform : streets to roads ids
print("Inverse transform : streets to roads ids")
start = time.time()
links_inverse_transform = {"STREET_ID":[],'GROUP_ID':[],"ROADS_IDS":[],'ROADS_SHARES':[]}
streets_list = list(Street_links["STREET_ID"])
links_inverse_transform["STREET_ID"] = streets_list
n_streets = len(streets_list)
links_inverse_transform['GROUP_ID'] = [[] for i in repeat(None, n_streets)]
links_inverse_transform['ROADS_IDS'] = [[] for i in repeat(None, n_streets)]
links_inverse_transform['ROADS_SHARES'] = [[] for j in repeat(None, n_streets)]
n_roads = len(links_transform["ROAD_ID"])
for r in range(n_roads) :
    road_r = links_transform["ROAD_ID"][r]
    matching_streets = links_transform["STREET_IDS"][r]
    shares = links_transform["SHARES"][r]
    group_g = links_transform["GROUP_ID"][r]
    for s, street_s in enumerate(matching_streets) :
        street_index = streets_list.index(street_s)
        links_inverse_transform["ROADS_IDS"][street_index].append(road_r)
        links_inverse_transform["ROADS_SHARES"][street_index].append(shares[s])
        links_inverse_transform['GROUP_ID'][street_index].append(group_g)
duration = time.gmtime(time.time()-start)
print("Inverse link transform computed... Time:", duration)

# 2) Distribution : (a) répartition des routes vers les rues et (b) des rues vers cellules
Street_links['GEOMETRY'] = Street_links['GEOMETRY'].apply(wkt.loads)
streets = gpd.GeoDataFrame(Street_links,crs='epsg:2154', geometry = 'GEOMETRY')
street_ids = {'STREET_ID':[], 'CELLS':[], 'DISTRIBUTION':[]}
street_ids['STREET_ID'] = streets_list
street_ids['CELLS'] = [[] for i in repeat(None, n_streets)]
start = time.time()
for poly_num, poly in enumerate(grid['polygon']):
    # Get the cell information
    x = grid['we_id'][poly_num]
    y = grid['sn_id'][poly_num]
    ## one file per cell with all roads included in it
    prefix = 'we_id_from_' + str(x) + '_sn_id_from_' + str(y)
    foutname = outdir + 'sub_' + prefix + '.dat'
    # Select roads intersecting the cell
    subdata = streets[streets.geometry.intersects(poly)]
    print('Road data croped in cell...')
    # For each road in the cell, get roads data and save cell number in the network file
    streets_in_cell = {'STREET_ID': list(subdata.STREET_ID), 'XA': list(subdata.XA), 'XB': list(subdata.XB),
                       'YA': list(subdata.YA), 'YB': list(subdata.YB), 'TOTAL_WIDTH':list(subdata.TOTAL_WIDTH),
                       'GEOMETRY':list(subdata.GEOMETRY), 'TUNNEL':list(subdata.TUNNEL)}
    for street_id in streets_in_cell['STREET_ID'] :
        street_index = Street_links[Street_links.STREET_ID == street_id].index[0] # street index in network file
        cells_list = street_ids['CELLS'][street_index] # get existing list of cells for this road
        cells_list.append(poly_num)
        street_ids['CELLS'][street_index] = cells_list
    # Save it
    print('Now saving....')
    pickle.dump(streets_in_cell, open(foutname, 'wb'))
    #Print information
    if poly_num % 100 == 0 :
        end = time.gmtime(time.time()-start)
        print("Distribution: done with cell ", poly_num, " in ", end, " for 100 cells" )
        start = time.time()
print("Streets in cells...")

for i, road_cells in enumerate(street_ids['CELLS']):
    n_cells = len(road_cells)
    distribution_list = [0] * n_cells
    street_ids['DISTRIBUTION'].append(distribution_list)
print('Distribution lists prepared...')

# for each street, calculate the percentage attributed to each cell it intersects
start = time.time()
for poly_num, poly in enumerate(grid['polygon']):
    # Select data in the cell
    x = grid['we_id'][poly_num]
    y = grid['sn_id'][poly_num]
    prefix = 'we_id_from_' + str(x) + '_sn_id_from_' + str(y)
    cell_streets_f = outdir + 'sub_' + prefix + '.dat'
    cell_streets_data = pickle.load(open(cell_streets_f,'rb'))
    cell_streets_df = pd.DataFrame(cell_streets_data)
    print('Data loaded...')
    # Calculate the proportion of the street in the cell for each street.
    for street_id in cell_streets_df['STREET_ID'] :
        street_index = Street_links[Street_links.STREET_ID == street_id].index[0]
        cells_list = street_ids['CELLS'][street_index]
        if len(cells_list) == 1 :
            street_ids['DISTRIBUTION'][street_index][0] = 1
        else :
            geom = Street_links.loc[street_index,'GEOMETRY']
            line_length = Street_links.loc[street_index,'LENGTH']
            intersection = poly.intersection(geom)
            intersection = sp_ops.transform(my_transformer.transform, intersection)
            intersection_length = intersection.length
            distrib = intersection_length / line_length
            i = cells_list.index(poly_num)
            street_ids['DISTRIBUTION'][street_index][i] = distrib
    end = time.gmtime(time.time()-start)
    print("Cell number ", poly_num, " - time ", end)

# Reassemble information
network_cells = pd.DataFrame(street_ids)
final_streets_df = pd.merge(left=Street_links, right=network_cells, how='left', on=['STREET_ID'])

# Save data
links_inverse_transform_df = pd.DataFrame.from_dict(links_inverse_transform)
links_inverse_transform_df.to_csv(links_inverse_transform_f, index=False)
final_streets_df.to_csv(street_links_attributes_f, index=False)
network_cells.to_csv(street_cells_distrib_f, index=False)

# next output dir
outdir = temp_dir + 'street_graph/attributes/subdomains_streets_height_width/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

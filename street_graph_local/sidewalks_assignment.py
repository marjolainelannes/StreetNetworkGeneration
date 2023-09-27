##################################################################################
# Study: Street network generation
# Purpose: Match sidewalks with streets in Paris
# Author: Marjolaine Lannes
# Creation date: December 12, 2022
##################################################################################
import sys
sys.path.append("../../include")
import pickle
import numpy as np
import pandas as pd
from shapely.geometry import LineString
from include.create_poly_around_road import create_1_poly_around_road

# Input and output directories
path = "C:/Users/lannesadm/PycharmProjects/Road_network_topology/"
sw_f = path + 'temp/load_sidewalks/Paris_sidewalks.dat'
streets_f = path + "temp/street_graph/graph_transform/street_links.dat"
output_f = path + "temp/street_graph/attributes/Paris_sidewalks.dat"

# Load data : streets & sidewalks
streets = pickle.load(open(streets_f, "rb"))
# streets = pd.read_csv(streets_f)
#streets['GEOMETRY'] = streets['GEOMETRY'].apply(wkt)
sidewalks_f = open(sw_f, 'rb')
sidewalks = pickle.load(sidewalks_f)
sidewalks_f.close()

# Create a dictionary of information for each street
print("initialize sw")
street_sidewalks = {"STREET_ID":[], 'SIDEWALKS': [], 'sw_widths': [], 'intersections': [], 'MEAN_SW_WIDTH': []}
for street_id in streets['STREET_ID']:
    street_sidewalks["STREET_ID"].append(street_id)
    street_sidewalks["SIDEWALKS"].append([])
    street_sidewalks["sw_widths"].append([])
    street_sidewalks["intersections"].append([])
    street_sidewalks["MEAN_SW_WIDTH"].append(None)
print("number of sidewalks:", len(sidewalks['LINE']))

# For each sidewalk, add its characteristics to the output data for the streets (in Paris) intersecting it.
print("Loop sw")
for i, line in enumerate(sidewalks['LINE']):
    if i%1000 == 0:
        print(i)
    # Get sidewalk id, summits and width (corrected by median width if error)
    sw_id = sidewalks['ID'][i]
    x1, x2 = line.xy[0]
    y1, y2 = line.xy[1]
    sw_width = sidewalks['WIDTH'][i]
    if sw_width == 0.0 or np.isnan(sw_width) or sw_width is None:
        sw_width = 3.1  # median sidewalk width in Paris
    poly = create_1_poly_around_road(x1, y1, x2, y2, sw_width/2 + 3 )
    found = False
    # Trial 1: if any street intersects the sidewalk polygon, add sidewalk characteristics to the road output data
    for j, l2 in enumerate(streets['GEOMETRY']):
        if poly.intersects(l2):
            street_id = streets.loc[j, 'STREET_ID']
            street_index = street_sidewalks['STREET_ID'].index(street_id)
            dx = poly.intersection(l2).length
            street_sidewalks['SIDEWALKS'][street_index].append(sw_id)
            street_sidewalks['sw_widths'][street_index].append(sw_width)
            street_sidewalks['intersections'][street_index].append(dx)
            found = True
    # Trial 2: try a larger polygon
    if not found:
        poly = create_1_poly_around_road(x1, y1, x2, y2, sw_width/2 + 6)
        for j, l2 in enumerate(streets['GEOMETRY']):
            if poly.intersects(l2):
                street_id = streets.loc[j, 'STREET_ID']
                street_index = street_sidewalks['STREET_ID'].index(street_id)
                street_sidewalks['SIDEWALKS'][street_index].append(sw_id)
                dx = poly.intersection(l2).length
                street_sidewalks['sw_widths'][street_index].append(sw_width)
                street_sidewalks['intersections'][street_index].append(dx)
                found = True
    # Trial 3: try a larger polygon
    if not found:
        poly = create_1_poly_around_road(x1, y1, x2, y2, sw_width/2 + 10)
        for j, l2 in enumerate(streets['GEOMETRY']):
            if poly.intersects(l2):
                street_id = streets.loc[j, 'STREET_ID']
                street_index = street_sidewalks['STREET_ID'].index(street_id)
                street_sidewalks['SIDEWALKS'][street_index].append(sw_id)
                dx = poly.intersection(l2).length
                street_sidewalks['sw_widths'][street_index].append(sw_width)
                street_sidewalks['intersections'][street_index].append(dx)
                found = True

# For each road intersected by sidewalks: averaging sidewalk information
print("Loop streets")
cnt = 0 ## initialize counter (of the number of streets with sidewalk data)
n_streets = len(list(street_sidewalks["STREET_ID"]))
for k in range(n_streets): ## for each street, calculate its mean sidewalk width
    if len(street_sidewalks['SIDEWALKS'][k]) > 0:
    # if there is at least 1 sidewalk in the street perimeter, calculated the weighted sw width
        cnt += 1
        tot = np.sum(street_sidewalks['intersections'][k])
        if tot > 0:
            aver = 0
            ws = street_sidewalks['intersections'][k] / tot
            for i, lg in enumerate(street_sidewalks['sw_widths'][k]):
                aver += lg * ws[i]
            street_sidewalks['MEAN_SW_WIDTH'][k] = aver

# Keep only mean sidewalk width and the list of sidewalk ids for each street
outdata = {'STREET_ID': street_sidewalks['STREET_ID'], 'SIDEWALKS': street_sidewalks['SIDEWALKS'],
           'MEAN_SW_WIDTH': street_sidewalks['MEAN_SW_WIDTH']}

# Save sidewalks data for this cell
fout = open(output_f,'wb')
pickle.dump(outdata, fout)
fout.close()

##################################################################################
# Study: Exposure modelling framework
# Purpose: Save Paris sidewalks sidewalks_data as shapefile
# Author: Myrto Valari
# Creation date: December 12, 2022
##################################################################################
from shapely.geometry import LineString
import pandas as pd
import geopandas as gpd

# Input and output files
path = "C:/Users/marjolaine.lannes/PycharmProjects/Road_network_topology/"
sidewalks_inputf = path + 'data/Paris_sidewalks/trottoirs.txt'
output_file = path + 'data/Paris_sidewalks/sidewalks_lamb93.shp'

# Initialize
myids = []
widths = []
lines = []
herids = []
herids1 = []

# Get sidewalks sidewalks_data in a Dataframe
with open(sidewalks_inputf) as f:
    for i, data in enumerate(f):
        if i > 0:
            myid = int(data.split(';')[0])
            herid = int(float(data.split(';')[1]))
            herid1 = int(float(data.split(';')[2]))
            width = float(data.split(';')[3])
            xa, xb, ya, yb = [float(x) for x in data.split(';')[-4:]]
            l = LineString(([xa, ya], [xb, yb]))
            myids.append(myid)
            herids.append(herid)
            herids1.append(herid1)
            widths.append(width)
            lines.append(l)
df = pd.DataFrame({'ID': myids, 'HERID': herids, 'HERID1': herids1, 'WIDTH': widths})

# Save it as shapefile in Lambert 93
gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=lines)
gdflamb93 = gdf.to_crs("EPSG:2154")
gdflamb93.to_file(output_file)

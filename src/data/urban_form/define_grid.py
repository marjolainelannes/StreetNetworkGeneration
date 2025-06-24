###################################################################################
# Study: Exposure modelling framework
# Purpose: Create a 0.02°x0.02° grid
# Note: The grid is defined in latlon projection = "epsg:4326"
###################################################################################
# documentation: https://fiona.readthedocs.io/en/latest/fiona.html#module-fiona
import pandas as pd
from shapely.geometry import Polygon
import pickle
import matplotlib.pyplot as plt
from pyproj import Proj, transform

path = "../../../"
outputdir = path + "data/grid_data/"

################################################################
# 1. Create a grid

def create_grid_from_coordinates(minx: float, dx: float, Nx: int, miny: float, dy: float, Ny: int, path: str):
   import fiona
   from shapely.geometry import mapping, LineString, MultiLineString
   # Set output files
   out_shp = path + 'grid.shp'
   # Set grid definition: starting from cell center
   maxx = minx + dx * (Nx - 1)
   maxy = miny + dy * (Ny - 1)
   # Create grid
   lines = []
   for x in range(Nx + 1):
      lines.append(
         LineString([(minx + (x - 1 / 2) * dx, miny - 1 / 2 * dy), (minx + (x - 1 / 2) * dx, maxy + 1 / 2 * dy)]))
   for y in range(Ny + 1):
      lines.append(
         LineString([(minx - 1 / 2 * dx, miny + (y - 1 / 2) * dy), (maxx + 1 / 2 * dx, miny + (y - 1 / 2) * dy)]))
   grid = MultiLineString(lines)
   # Save grid as shapefile
   schema = {
      "geometry": "MultiLineString",
      "properties": {"id": "int"}
   }
   with fiona.open(out_shp, 'w', driver="ESRI Shapefile", schema=schema, crs="epsg:4326") as grid_shp:
      grid_shp.write({'geometry': mapping(grid), "properties": {"id": 0}})

create_grid_from_coordinates(minx=1.35, dx=0.02, Nx=110, miny=48, dy=0.02, Ny=75, path=outputdir)

###############################################################
# 2. Create polygons for cells

def create_cell_polys(minx, dx, Nx, miny, dy, Ny, buffer, path) :
   fnameout_no_buffer='polygons_cells_IDF_no_buffer_lambert93.dat'
   fnameout_with_buffer = 'polygons_cells_IDF_with_buffer_'+str(buffer)+'m_lambert93.dat'
   fnameout_lonlat='polygons_cells_IDF_no_buffer_lonlat.dat'
   fnameout_csv = 'polygons_cells_IDF_no_buffer_lonlat.csv'
   lamb93=Proj(init='epsg:2154')
   lonlat=Proj(init='epsg:4326')

   # for each cell, get lat, lon ids and cell polygon
   polys_latlon = {'we_id': [], 'sn_id': [], 'polygon': []}
   polys_lamb93 = {'we_id': [], 'sn_id': [], 'polygon': []}
   polys_lamb93_with_buffer = {'we_id': [], 'sn_id': [], 'polygon': []}
   for y in range (0,Ny):
      for x in range (0,Nx):
         #### Latlon
         # Define summits for each cell in latlon
         lon_bl = minx + (x - 1/2)*dx
         lat_bl = miny + (y - 1/2)*dy
         lon_br = minx + (x + 1/2)*dx
         lat_br = miny + (y - 1/2)*dy
         lon_tr = minx + (x + 1 / 2) * dx
         lat_tr = miny + (y + 1 / 2) * dy
         lon_tl = minx + (x - 1 / 2) * dx
         lat_tl = miny + (y + 1 / 2) * dy
         # lon_bl = minx + dx * x
         # lat_bl = miny + dy * y
         # lon_br = minx + dx * (x+1)
         # lat_br = miny + dy * y
         # lon_tr = minx + dx * (x+1)
         # lat_tr = miny + dy * (y+1)
         # lon_tl = minx + dx * x
         # lat_tl = miny + dy * (y+1)
         # save (x,y) ids and the cell in polys_latlon dictionnary
         poly = Polygon([(lon_bl, lat_bl), (lon_br, lat_br), (lon_tr, lat_tr), (lon_tl, lat_tl), (lon_bl, lat_bl)])
         polys_latlon['we_id'].append(x)
         polys_latlon['sn_id'].append(y)
         polys_latlon['polygon'].append(poly)
         #### Lambert 93
         # Define summits for each cell in lambert93
         x_bl, y_bl = transform(lonlat, lamb93, lon_bl, lat_bl)
         x_br, y_br = transform(lonlat, lamb93, lon_br, lat_br)
         x_tr, y_tr = transform(lonlat, lamb93, lon_tr, lat_tr)
         x_tl, y_tl = transform(lonlat, lamb93, lon_tl, lat_tl)
         # Save the grid file (without buffer)
         poly_no_buffer = Polygon([(x_bl, y_bl), (x_br, y_br), (x_tr, y_tr), (x_tl, y_tl)])
         polys_lamb93['we_id'].append(x)
         polys_lamb93['sn_id'].append(y)
         polys_lamb93['polygon'].append(poly_no_buffer)
         # Save the grid file (with buffers)
         poly_with_buffer = Polygon([(x_bl-buffer, y_bl-buffer), (x_br+buffer, y_br-buffer), (x_tr+buffer, y_tr+buffer), (x_tl-buffer, y_tl+buffer)])
         polys_lamb93_with_buffer['we_id'].append(x)
         polys_lamb93_with_buffer['sn_id'].append(y)
         polys_lamb93_with_buffer['polygon'].append(poly_with_buffer)
         # Plot the grid
         plt.plot(*poly_with_buffer.exterior.xy, 'k')
   pickle.dump(polys_latlon, open(path + fnameout_lonlat, 'wb'))
   pickle.dump(polys_lamb93, open(path + fnameout_no_buffer, 'wb'))
   pickle.dump(polys_lamb93_with_buffer, open(path + fnameout_with_buffer, 'wb'))
   polys_df = pd.DataFrame.from_dict(polys_latlon)
   polys_df.to_csv(fnameout_csv, index=False)
   plt.show()

create_cell_polys(minx=1.35, dx=0.02, Nx=110, miny=48, dy=0.02, Ny=75, buffer=100, path=outputdir)

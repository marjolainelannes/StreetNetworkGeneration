# Street network generation model



This project builds a street network with street length, width and mean building's height based on OpenStreetMap road network and a building database.  
Its original use is a street geometry preprocessor for street air quality models such as MUNICH (Kim et al. 2022). 

## Model description

![Flow_diagram](docs/Flow_diagram.png "Flow_diagram")

## Input data
  
Mandatory input files:
+ OpenStreetMap road network: here the xml output from eqasim (based on pt2matsim) is read  
+ buildings: BDNB database  
  
Optional input file:  
+ OpenStreetMap detailed road network (with detailed geometry instead of a straight line): here the xml output from eqasim (based on pt2matsim) is read  

## Simulation

This section describes command lines to run the model. Some codes are parallelized with algo. Calculation times for Paris region (600,000 links) are given when it takes more than a few minutes.

1 - Load data describing the urban form in data/urban_form  
* python define_grid.py  
* python adjacent_cells.py 
* python load_BDNB.py
* python makeSUBbuildings.py (1h14)  
  
2 - Load OpenStreetMap data in data/load_OSM  
* python read_OSM.py
* python road_graph.py (16h40)  
* python road_links_neighborhood (1h30)  
* python makeSUBroads.py (11h)  
* python load_waterways.py  
  
3 - Build the street graph  
* python groups.py (3d15h)  
* python graph_transform.py (1d16h)  
* python streets_distribution.py (10h)  
* create the output folder, then run the following (1h):  
./algo start --computer-file=nodes_urbino.dat --argument-file=Table_of_cells_urbino.dat --log=log_urbino.txt python run_street_edges_height.py  
* python attributes_transform.py (1h30)

If you already have an appropriate street network and would like to add width and height attributes to street links, you shall use the following code: src/street_graph/street_edges_height.py  

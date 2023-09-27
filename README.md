# Street_network_generation_model

This project constructs street network topology from eqasim output.
It is meant to be used as a street geometry preprocessor for POLAIR/MUNICH.

## 1 - load BDTOPO data in data/load_BDTOPO
a) python makeSUBbuildings.py
b) ./algo start --computer-file=nodes_byblos.dat --argument-file=Table_of_cells_byblos.dat --log=log_SUBroads_byblos.txt python run_makeSUBroads.py

## 2 - load OpenStreetMap data in data/load_OSM
a) python read_OSM.py
b) python road_graph.py (10h)
c) ./algo start --computer-file=nodes_luni.dat --argument-file=Table_of_cells_luni.dat --log=log_SUBroads_luni.txt python run_makeSUBroads.py (1h)
d) python roads_in_cells.py

## 3 - load sidewalks in data/load_sidewalks
a) python save_Paris_sidewalks_as_shp.py
b) ./algo start --computer-file=nodes_carthage.dat --argument-file=Table_of_cells_carthage.dat --log=log_SUBsidewalks_carthage.txt python run_makeSUBsidewalks.py

## 4 - roads graph
a) road_width :
./algo start --computer-file=nodes_luni.dat --argument-file=Table_of_cells_luni.dat --log=log_luni.txt python run_road_width.py
b) road_graph_synthesis

## 5 - streets graph
a) groups
b) graph_transform
c) sidewalks_assignment
d) street_width_distribution
e) street_edges_heights
f) attributes_transform
g) roads_to_streets

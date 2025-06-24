##################################################################################
# Study: Street network generation
# Purpose: to save buildings information from BDNB
# Note: Get buildings information (id, height, geometry) from BDNB for Ile-de-France region
##################################################################################
import os
import pandas as pd
import pickle

# Config
path = "../../../"
input_dir = path + "data/"
cache_dir  = path + "temp/"

# Input / output files
buildings_dir = input_dir + 'BDNB/'
foutname = cache_dir + 'data/urban_form/buildings_BDNB.dat'
foutname_csv = cache_dir + 'data/urban_form/buildings_BDNB.csv'
f = open(foutname, "wb")

# Load data
buildings = pd.DataFrame(columns=['ID','GEOMETRY','HEIGHT'])
for buildings_f in os.listdir(buildings_dir) :
    ## Load buildings for each department
    print(buildings_f)
    buildings_department = pd.read_csv(buildings_dir + buildings_f)
    ## Filter fictive buildings
    if "fictive_geom_cstr" in buildings_department.columns :
        buildings_department_filtered = buildings_department[(buildings_department.fictive_geom_cstr == 0)]
    else :
        buildings_department_filtered = buildings_department
    print(buildings_department_filtered)
    ## Get height information
    buildings_department_filtered = buildings_department_filtered[['WKT','batiment_construction_id','hauteur']]
    buildings_department_cleaned = buildings_department_filtered.rename(columns={'batiment_construction_id':'ID','WKT':'GEOMETRY','hauteur':'HEIGHT'})
    print(buildings_department_cleaned)
    buildings = pd.concat([buildings,buildings_department_cleaned])
    print(buildings)

# Save outputs
print(buildings)
buildings.to_csv(foutname_csv, index=False)
pickle.dump(buildings, f)
f.close()


###################################################################################
# Study: Street network generation
# Purpose: Get adjacent cells
###################################################################################
import pickle

# Config
path = "../../../"
input_dir = path + "data/"
Nx=110
Ny=75

# Get a list of adjacent cells for each cell
N_cells = Nx * Ny
adjacent_cells = {}
for C in range(N_cells):
    # We write C = y * Nx + x
    y = int(C/Nx)
    x = C % Nx
    # do we have to add right/left cells ?
    x_iterate = []
    if x > 0 :
        x_iterate.append(x-1)
    x_iterate.append(x)
    if x < Nx-1 :
        x_iterate.append(x+1)
    # add upper cells
    adjacent_cell_C = []
    if (y > 0)  :
        for k in x_iterate:
            adjacent_cell_C.append((y-1) * Nx + k)
    # add cells on the same longitude
    for k in x_iterate:
        adjacent_cell_C.append(y * Nx + k)
    # add bottom cells
    if (y < Ny-1):
        for k in x_iterate:
            adjacent_cell_C.append((y + 1) * Nx + k)
    # add cell C adjacency list to dict
    adjacent_cells[C] = adjacent_cell_C

# Save it
f_out = input_dir + "grid_data/adjacent_cells.dat"
pickle.dump(adjacent_cells, open(f_out, "wb"))

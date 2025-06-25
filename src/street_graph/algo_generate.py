#######################################################################
# Note: generate lists of cells and nodes for each computer to run the code in parallel
#######################################################################
import pandas as pd

## Parameters
N_cells = 8250
N_processors = 8

## Table of cells
output_1 = "Table_of_cells_1.dat"
output_2 = "Table_of_cells_2.dat"
n = N_cells/2
list_1 = []
list_2 = []
for i in range (0,n):
    list_1.append(i)
    list_2.append(n+i)
df_1 = pd.DataFrame({'NCELL': list_1})
df_1.to_csv(output_1, index=False)
df_2 = pd.DataFrame({'NCELL': list_2})
df_2.to_csv(output_2, index=False)

## Nodes
computers = ['computer_1','computer_2']
for c in computers :
    output = "nodes_" + c + ".dat"
    content = [str(N_processors) + "/" + c]
    df = pd.DataFrame(content)
    df.to_csv(output, index=False, header=False)

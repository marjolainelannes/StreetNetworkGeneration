import pandas as pd

### Tables of cells

# Output files
output_cathage = "Table_of_cells_carthage.dat"
output_halong = "Table_of_cells_halong.dat"
# Add cells
list_c = []
list_h = []
for i in range (0,4125):
    list_c.append(i)
for i in range (4125,8250):
    list_h.append(i)
# Save outputs
df_c = pd.DataFrame({'NCELL': list_c})
df_c.to_csv(output_cathage, index=False)
df_h = pd.DataFrame({'NCELL': list_h})
df_h.to_csv(output_halong, index=False)

## Nodes
machines = ["carthage", "halong"]
for m in machines :
    output = "nodes_" + m + ".dat"
    contenu = ["3/" + m]
    df = pd.DataFrame(contenu)
    df.to_csv(output, index=False, header=False)
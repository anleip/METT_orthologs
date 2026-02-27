# %%
import numpy as np
import pandas as pd

# %%
# Filter the table according to the wet lab team
all_genes = pd.read_csv("/nfs/research/jlees/anlei/analysis/METT_ortholog/essential_genes_uniformis_vulgatus_table_only.csv")
print(all_genes)
# %%
essential_genes = all_genes[(all_genes['media']=="solid")&(all_genes['category']=="conserved_ES")]
print(essential_genes)
# %%
bu_essential_genes = essential_genes[essential_genes['species']=="uniformis"]
print(bu_essential_genes)
# %%

# %%
from proscope.protein import Protein
a = Protein('MYC')
# %%
a.low_or_high_plddt_region
# %%
ax = a.plot_plddt(filename='../data/myc_protein_plot.png')
# save plot to data folder
# %%
a.low_or_high_plddt_region_sequence
# %%
a.domains
# %%

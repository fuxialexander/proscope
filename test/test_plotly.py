# %%
from proscope.protein import Protein
a = Protein('PAX5', use_es=False, window_size=10)
#%%
a.plotly_plddt()
# %%

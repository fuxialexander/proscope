# %%
from proscope.protein import Protein
import pandas as pd 
import numpy as np

a = Protein('TP53', window_size=10, use_es = True, af2_folder='./', esm_folder='~/Projects/esm/examples/variant-prediction/output_cancer_genes')
#%%
afmis = pd.read_feather("./AlphaMissense_aa_substitutions.feather")
#%%
af_g = afmis.query('uniprot_id==@a.uniprot_id')
af_g['pos'] = af_g['protein_variant'].str[1:-1].astype(int)
af_g['aa'] = af_g['protein_variant'].str[-1]
af_g_wide = af_g.pivot(index='pos', columns='aa', values='am_pathogenicity')
#%%
from proscope.protein import smooth
to_compare = {'AFMissense': ('lime', smooth(af_g_wide.mean(1).values/af_g_wide.mean(1).values.max(), 10)), 
              # 'ESM1b': ('lightpink', a.esm),
              # 'pLDDT_original': ('blue', a.plddt/100),
              'ES': ('green', a.es_raw)
              }
# %%
a.plotly_plddt(to_compare=to_compare, domains_to_show=['DNA-binding region'])


# %%z
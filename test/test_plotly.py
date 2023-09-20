# %%
from proscope.protein import Protein
import pandas as pd 
import numpy as np

a = Protein('PAX5', window_size=10, use_es = True, af2_folder='./', esm_folder='~/Projects/proscope/data/content/ALL_hum_isoforms_ESM1b_LLR')
#%%
esm1v = pd.read_feather("./all_gene_dfs.feather")
# %%
PAX5 = esm1v.query('Gene == "PAX5"')
# %%
# PAX5 to wide df
PAX5['pos'] = PAX5['Mutation'].str[1:-1].astype(int)
PAX5['aa'] = PAX5['Mutation'].str[-1]
PAX5_wide = PAX5.pivot(index='pos', columns='aa', values='esm1v_t33_650M_UR90S_1')
# %%
afmis = pd.read_feather("./AlphaMissense_aa_substitutions.feather")
PAX5_af = afmis.query('uniprot_id==@a.uniprot_id')
PAX5_af['pos'] = PAX5_af['protein_variant'].str[1:-1].astype(int)
PAX5_af['aa'] = PAX5_af['protein_variant'].str[-1]
PAX5_af_wide = PAX5_af.pivot(index='pos', columns='aa', values='am_pathogenicity')
#%%
to_compare = {'AFMissense': ('lime', PAX5_af_wide.mean(1).values), 
              'ESM1v': ('purple', PAX5_wide.mean(1).values/PAX5_wide.mean(1).values.min()),
              'ESM1b': ('lightpink', a.esm),
            #   'pLDDT_original': ('blue', a.plddt/100),
              'ES': ('green', a.es)
              }
# %%
a.plotly_plddt(to_compare=to_compare)

# %%
a.plot_plddt_manuscript()
# %%

# %%
from proscope.af2 import AFPairseg

# %%
a = AFPairseg(
    "/home/xf2217/Projects/proscope/test/test_causal_db/structures/causal/MECP2_TFAP2A",
    "/home/xf2217/Projects/proscope/test/test_causal_db/sequences/causal_edges_db/MECP2_TFAP2A")
# %%
a.plotly_plddt_gene1()
# %%
a.plotly_plddt_gene2()
# %%
a.plot_plddt_gene1()
#%%
a.plot_plddt_gene2()
# %%
a.protein1.plotly_plddt()
# %%
a.protein2.plotly_plddt()
# %%
a.plot_score_heatmap()
# %%
from proscope.viewer import view_pdb_html
view_pdb_html(a.pairs_data['TAF1_3_ZFX_0'].pdb)

# %%
from proscope.af2 import read_pdb
read_pdb('/home/xf2217/Projects/proscope/test/test_causal_db/structures/causal/MECP2_TFAP2A/MECP2_1_TFAP2A_1/MECP2_1.TFAP2A_1_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_042.pdb')
# %%

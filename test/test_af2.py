# %%
from proscope.af2 import AFPairseg

# %%
a = AFPairseg("/home/xf2217/Projects/proscope/test/test_causal_db/structures/causal/TAF1_ZFX",
          "/home/xf2217/Projects/proscope/test/test_causal_db/sequences/causal_edges_db/TAF1_ZFX")
# %%
a.plot_plddt_gene1()
# %%
a.plot_plddt_gene2()
# %%
a.protein1.plot_plddt()
# %%
a.protein2.plot_plddt()
# %%
a.plot_score_heatmap()
# %%

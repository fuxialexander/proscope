# %%
# use nglview to visualize the structure
import nglview as nv
# %%
# homodimer structure align with the above
homodimer = 'test_causal_db/structures/homodimer/PRDM1/PRDM1_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb'
dimer = '/home/xf2217/Projects/proscope/test/test_causal_db/structures/causal/PRDM1_SMAD2/PRDM1_4_SMAD2_3/PRDM1_4.SMAD2_3_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_042.pdb'
# %%
# open both homodimer and dimer using nv.FileStructure and view in the same widget
view = nv.NGLWidget()
view.add_component(nv.FileStructure(homodimer))
view.add_component(nv.FileStructure(dimer))
# color the homodimer in red and dimer in blue
view.update_cartoon(color='red', component=0)
view.update_cartoon(color='blue', component=1)
# align the dimer to the homodimer 
view.align(component=1, ref_component=0)
view
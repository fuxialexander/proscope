# %%
import sys
# sys.path.append('/home/xf2217/Projects/proscope/proscope')
from proscope.protein import Protein, generate_pair_sequence
import os
#%%
generate_pair_sequence('PAX5', 'RARA', 'test_causal_db/sequences/causal_edges_db/PAX5_RARA')
#%%
generate_pair_sequence('TAL1', 'SIN3A', 'test_causal_db/sequences/other/')
#%%
import pandas as pd
df = pd.read_csv('../data/content/ALL_hum_isoforms_ESM1b_LLR/O95274_LLR.csv', index_col=0)
# melt to long format, column, row, value
df = df.reset_index().melt(id_vars='index')
df['variant'] = df['variable'].str.replace(' ', '') + df['index'].astype(str)
df['pos'] = df['variable'].apply(lambda x: int(x.split(' ')[1]))
df = df.rename({'value': 'esm'}, axis=1)
df = df[['variant', 'pos', 'esm']]
#%%
a = Protein('PAX5')
#%%
fig, ax = a.plot_plddt(show_domain=True, show_low_plddt=False)
# add legent
ax.set_xlim(0, len(a.plddt))
ax.set_xticks(range(0, len(a.plddt), 500))
# add lines at 300, 400, 500 with height 0.1, 0.3, 0.5
ax.axvline(30, ymax=0.1, color='grey', linestyle='-')
ax.axvline(40, ymax=0.3, color='grey', linestyle='-')
ax.axvline(50, ymax=0.5, color='grey', linestyle='-')
#%%
b = Protein('STAT3')
b.plot_plddt()
# %%
a.low_or_high_plddt_region
# %%
ax = a.plot_plddt(filename='../data/rela_protein_plot.png')
# save plot to data folder
# %%
a.low_or_high_plddt_region_sequence
# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
def plot_plddt(obj, to_compare=None, filename=None):
    plt.figure(figsize=(12, 3))
    plt.plot(obj.plddt)
    if to_compare is not None:
        plt.plot(to_compare)
    # highlight low plddt region
    for region in obj.low_plddt_region:
        plt.axvspan(region[0], region[1], ymax=1, ymin=0, color='grey', alpha=0.1)
    
    # highlight domain, color by feature_type
    cmap = plt.get_cmap('tab20').colors
    # map feature_type to color
    feature_type_to_color = {}
    obj.domains = obj.domains.query('(feature_type=="domain") or (feature_type=="region of interest")')
    # for i, t in enumerate(obj.domains.feature_type.unique()):
    #     feature_type_to_color[t] = cmap[i+5]
    feature_type_to_color['domain'] = cmap[5]
    feature_type_to_color['region of interest'] = cmap[6]
    

    y_span = 0.1#0.8/len(obj.domains.feature_type.unique())
    for i, domain in obj.domains.iterrows():
        idx = np.where(obj.domains.feature_type.unique()==domain.feature_type)[0][0]
        plt.axvspan(domain.feature_begin, domain.feature_end, ymax=idx * y_span + y_span, ymin=idx * y_span 
            , color=feature_type_to_color[domain.feature_type], alpha=0.2)
    # add legend of domain color
    legend_elements = []
    for i in obj.domains.feature_type.unique():
        legend_elements.append(Patch(facecolor=feature_type_to_color[i], label=i))
        #reverse the order of legend
    legend_elements = legend_elements[::-1]
    # add "low plddt region" to legend
    legend_elements.append(Patch(facecolor='grey', label='low pLDDT region'))
    # add number index of low or high plddt region on top of the plot
    # for i, region in enumerate(obj.low_or_high_plddt_region):
    #     plt.text(region[0], 0.9, f"{i}", fontsize=12)
    # legend outside the plot
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.ylabel(f"{obj.gene_name} pLDDT")
    plt.xlabel("Residue")
    # set xlim to the length of protein
    plt.xlim(0, len(obj.plddt))
    # plt.ylabel("pLDDT")
    # plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

plot_plddt(a)
# %%
b = Protein("NKX2-1")
# %%
a.plot_plddt(
    
)
# %%
b.plot_plddt()

#%%
def check_sequence_similarity(a,b):
    """check sequence similarity using pairwise2"""
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    alignments = pairwise2.align.globalxx(a, b)
    print(format_alignment(*alignments[0]))
    return alignments[0][2]*2/(len(a)+len(b))

x = str(a.low_or_high_plddt_region_sequence[1].seq)
y = str(b.low_or_high_plddt_region_sequence[1].seq)
check_sequence_similarity(x,y)
# %%
str(a.low_or_high_plddt_region_sequence[1].seq) + ':' + str(b.low_or_high_plddt_region_sequence[0].seq)
# %%
str(a.sequence) + ':' + str(a.sequence) + ':' + str(b.sequence)
# %%
c = Protein('EP300')
# %%
c.plot_plddt()
# %%
for i in [1,3,7,8,10]:
    with open(f"snai_rela_ep300_{i}.fasta", 'w') as f:
        print(f">snai_rela_ep300_{i}", file=f)
        print(str(a.low_or_high_plddt_region_sequence[1].seq) + ':' + str(b.low_or_high_plddt_region_sequence[0].seq) + ':' + str(c.low_or_high_plddt_region_sequence[i].seq), file=f)
# %%
df = Protein('TFAP2C')
# %%
b.plot_plddt()
# %%
c.low_or_high_plddt_region_sequence[1]
# %%
str(a.low_or_high_plddt_region_sequence[2].seq) + ':' + \
      str(b.sequence[0:310]) 
# %%
for i,s in enumerate(a.low_or_high_plddt_region_sequence):
    for j,s1 in enumerate(b.low_or_high_plddt_region_sequence):
        with open(f"zeb1_{i}_esr1_{j}.fasta", 'w') as f:
            print(f">zeb1_{i}_esr1_{j}", file=f)
            print(str(s.seq) + ':' + str(s1.seq),  file=f)
# %%
str(a.sequence) + ':' + \
      str(b.sequence)  + ':' + \
        str(b.sequence[653:])  + ':' + \
            str(b.sequence[0:327])  + ':' + \
        str(b.sequence[653:]) 
# %%

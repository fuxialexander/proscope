# %%
import sys
# sys.path.append('/home/xf2217/Projects/proscope/proscope')
from proscope.protein import Protein, generate_pair_sequence
import os
#%%
generate_pair_sequence('TFAP2A', 'ZFX', 'test_causal_db/sequences/causal_edges_db/TFAP2A_ZFX')
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

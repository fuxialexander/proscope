# %%
import multiprocessing as mp
import os

import pandas as pd
from tqdm import tqdm

from proscope.protein import Protein

pairs = pd.read_csv('causal_top_pair_per_celltype_df_non_redundant_pair_with_count.csv', sep=',').query('count>5')

os.makedirs('causal_edges_db', exist_ok=True)

def process_row(arg):
    k, row = arg
    try:
        protein_a = Protein(row['source'])
        protein_b = Protein(row['target'])
        low_or_high_plddt_region_sequence_a = protein_a.low_or_high_plddt_region_sequence
        low_or_high_plddt_region_sequence_b = protein_b.low_or_high_plddt_region_sequence
        for i, seq_a in enumerate(low_or_high_plddt_region_sequence_a):
            for j, seq_b in enumerate(low_or_high_plddt_region_sequence_b):
                os.makedirs(f"causal_edges_db/{protein_a.gene_name}_{protein_b.gene_name}", exist_ok=True)
                filename = f"causal_edges_db/{protein_a.gene_name}_{protein_b.gene_name}/{protein_a.gene_name}_{i}_{protein_b.gene_name}_{j}.fasta"
                with open(filename, 'w') as f:
                    f.write(f">{protein_a.gene_name}_{i}.{protein_b.gene_name}_{j}\n{str(seq_a.seq)}:{str(seq_b.seq)}\n")
    except:
        pass
    
def parallel_process():
    num_processes = mp.cpu_count()-2
    pool = mp.Pool(num_processes)
    results = pool.map(process_row, pairs.iterrows())
    pool.close()
    pool.join()

if __name__ == '__main__':
    parallel_process()

# %%

# %%
from glob import glob
# %%
a3m = glob('*a3m')
# %%
# for f in a3m, read second line and split by space and take second element
id_map_file = {}
id_map_folder = {}
id_map_root = {}
for a in a3m:
    with open(a, 'r') as f:
        f.readline()
        print(a)
        id = a.split('.')[0]
        pair = f.readline().split()[1]
        id_map_file[id] = pair
        id_map_folder[id] = pair.replace('.', '_')
        tf1, _, tf2, _ = pair.replace('.', '_').split('_')
        id_map_root[id] = f'{tf1}_{tf2}'

# %%
id_map_root
# %%
# reorganize every file in the 'output' folder into the following structure in 'output_reorganized'
# from output/original_basename which looks like f'{id}_.*'
# output_reorganized/id_map_root[id]/id_map_folder[id]/original_basename.replace(id, id_map_file[id])
# use pure python
import os
for f in glob('output/[0-9]*'):
    original_basename = f.split('/')[-1]
    id = original_basename.split('.')[0].split('_')[0]
    new_basename = original_basename.replace(id, id_map_file[id])
    root = id_map_root[id]
    folder = id_map_folder[id]
    print(f'{root}/{folder}/{new_basename}')
    os.makedirs(f'output_reorganized/{root}/{folder}', exist_ok=True)
    os.rename(f, f'output_reorganized/{root}/{folder}/{new_basename}')
    # copy config.json to every folder
    os.system(f'cp output/config.json output_reorganized/{root}/{folder}/config.json')
# %%

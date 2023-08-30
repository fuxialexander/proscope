import os
import gzip
import pandas as pd
import numpy as np
import pickle
from tqdm import tqdm
from Bio import SeqIO
import xmlschema

_seq = None
_genename_to_uniprot = None
_lddt = None
_schema = None


def get_seq():
    global _seq
    if _seq is None:
        initialize_seq()
    return _seq

def get_genename_to_uniprot():
    global _genename_to_uniprot
    if _genename_to_uniprot is None:
        initialize_genename_to_uniprot()
    return _genename_to_uniprot

def get_lddt():
    global _lddt
    if _lddt is None:
        initialize_lddt()
    return _lddt

def get_schema():
    global _schema
    if _schema is None:
        initialize_schema()
    return _schema

def initialize_schema():
    global _schema
    _schema = xmlschema.XMLSchema('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd')


def initialize_seq():
    global _seq
    _seq = {}
    with gzip.open(f'{os.path.dirname(__file__)}/uniprot_sprot.fasta.gz', 'rt') as f:
        for record in tqdm(SeqIO.parse(f, "fasta")):
            id = record.id.split('|')[1]
            _seq[id] = record.seq

def initialize_genename_to_uniprot():
    global _genename_to_uniprot
    _genename_to_uniprot = pd.read_csv(f'{os.path.dirname(__file__)}/uniprot_to_genename.txt', sep='\t').set_index('To').to_dict()['From']

def initialize_lddt():
    global _lddt
    if not os.path.exists(f'{os.path.dirname(__file__)}/9606.pLDDT.pickle'):
        _lddt = dict()
        with open(f'{os.path.dirname(__file__)}/9606.pLDDT.tdt', 'r') as f:
            for line in f:
                id, score = line.strip().split('\t')
                _lddt[id] = np.array(score.split(",")).astype(float)
        with open(f'{os.path.dirname(__file__)}/9606.pLDDT.pickle', 'wb') as f:
            pickle.dump(_lddt, f)
    else:
        with open(f'{os.path.dirname(__file__)}/9606.pLDDT.pickle', 'rb') as f:
            _lddt = pickle.load(f)

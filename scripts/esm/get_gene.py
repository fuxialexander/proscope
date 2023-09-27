import argparse
import csv
import gzip
import itertools
import json
import logging
import os
import pathlib
import warnings
from argparse import ArgumentParser
from collections import defaultdict
from typing import List, Tuple

import numpy as np
import pandas as pd
import torch
from Bio import BiopythonWarning, SeqIO
from Bio.PDB import is_aa
from Bio.SeqUtils import seq1
from esm import (Alphabet, FastaBatchedDataset, MSATransformer,
                 ProteinBertModel, pretrained)
from tqdm import tqdm


def create_parser():
    parser = argparse.ArgumentParser(
        description="Label a deep mutational scan with predictions from an ensemble of ESM-1v models."  # noqa
    )
    # fmt: off
    parser.add_argument(
        "--genes_file",
        type=str,
        help="Genes to predict")
    parser.add_argument(
        "--gene",
        type=str,
        help="Gene to predict")
    parser.add_argument("--out_dir", type=str, help="Output directory", default=".")
    parser.add_argument(
        "--model_location",
        type=str,
        help="PyTorch model file OR name of pretrained model to download (see README for models)",
        nargs="+",
    )
    parser.add_argument(
        "--sequence",
        type=str,
        help="Base sequence to which mutations were applied",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        help="fasta file location",
    )
    parser.add_argument(
        "--dms-input",
        type=pathlib.Path,
        help="CSV file containing the deep mutational scan",
    )
    parser.add_argument(
        "--mutation-col",
        type=str,
        default="Mutation",
        help="column in the deep mutational scan labeling the mutation as 'AiB'"
    )
    parser.add_argument(
        "--dms-output",
        type=pathlib.Path,
        help="Output file containing the deep mutational scan along with predictions",
    )
    parser.add_argument(
        "--offset-idx",
        type=int,
        default=0,
        help="Offset of the mutation positions in `--mutation-col`"
    )
    parser.add_argument(
        "--scoring-strategy",
        type=str,
        default="wt-marginals",
        choices=["wt-marginals", "pseudo-ppl", "masked-marginals"],
        help=""
    )
    parser.add_argument(
        "--msa-path",
        type=pathlib.Path,
        help="path to MSA (required for MSA Transformer)"
    )
    parser.add_argument(
        "--msa-samples",
        type=int,
        default=400,
        help="number of sequences to randomly sample from the MSA"
    )
    # fmt: on
    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
    return parser


def get_uniprot_accession(record):
    """

    Parameters
    ----------
    record

    Returns
    -------

    """
    parts = record.id.split('|')
    return parts[1]


def read_msa(filename: str, nseq: int) -> List[Tuple[str, str]]:
    """ Reads the first nseq sequences from an MSA file, automatically removes insertions."""

    msa = [
        (record.description, str(record.seq))
        for record in itertools.islice(SeqIO.parse(filename, "fasta"), nseq)
    ]
    msa = [(desc, seq.upper()) for desc, seq in msa]
    return msa

def label_row(row, sequence, token_probs, alphabet, offset_idx):
    breakpoint()
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence:" + sequence[idx] + str(idx) + "-" + wt

    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)

    # add 1 for BOS
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()


def compute_pppl(row, sequence, model, alphabet, offset_idx):
    wt, idx, mt = row[0], int(row[1:-1]) - offset_idx, row[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"

    # modify the sequence
    sequence = sequence[:idx] + mt + sequence[(idx + 1) :]

    # encode the sequence
    data = [
        ("protein1", sequence),
    ]

    batch_converter = alphabet.get_batch_converter()

    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)

    # compute probabilities at each position
    log_probs = []
    for i in range(1, len(sequence) - 1):
        batch_tokens_masked = batch_tokens.clone()
        batch_tokens_masked[0, i] = alphabet.mask_idx
        with torch.no_grad():
            token_probs = torch.log_softmax(model(batch_tokens_masked.cuda())["logits"], dim=-1)
        log_probs.append(token_probs[0, i, alphabet.get_idx(sequence[i])].item())  # vocab size
    return sum(log_probs)


def get_output(df, seq, offset, models, args):
    for m in models:
        model, batch_converter, alphabet, model_location = m
        # inference for each model

        if isinstance(model, MSATransformer):
            data = [read_msa(args.msa_path, args.msa_samples)]
            assert (
                args.scoring_strategy == "masked-marginals"
            ), "MSA Transformer only supports masked marginal strategy"

            batch_labels, batch_strs, batch_tokens = batch_converter(data)

            all_token_probs = []
            for i in tqdm(range(batch_tokens.size(2))):
                batch_tokens_masked = batch_tokens.clone()
                batch_tokens_masked[0, 0, i] = alphabet.mask_idx  # mask out first sequence
                with torch.no_grad():
                    token_probs = torch.log_softmax(
                        model(batch_tokens_masked.cuda())["logits"], dim=-1
                    )
                all_token_probs.append(token_probs[:, 0, i])  # vocab size
            token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
            df[model_location] = df.apply(
                lambda row: label_row(
                    row[args.mutation_col], seq, token_probs, alphabet, offset
                ),
                axis=1,
            )

        else:
            data = [
                ("protein1", seq),
            ]
            batch_labels, batch_strs, batch_tokens = batch_converter(data)

            if args.scoring_strategy == "wt-marginals":
                with torch.no_grad():
                    token_probs = torch.log_softmax(model(batch_tokens.cuda())["logits"], dim=-1)
                df[model_location] = df.apply(
                    lambda row: label_row(
                        row[args.mutation_col],
                        seq,
                        token_probs,
                        alphabet,
                        offset,
                    ),
                    axis=1,
                )
            elif args.scoring_strategy == "masked-marginals":
                all_token_probs = []
                for i in tqdm(range(batch_tokens.size(1))):
                    batch_tokens_masked = batch_tokens.clone()
                    batch_tokens_masked[0, i] = alphabet.mask_idx
                    with torch.no_grad():
                        token_probs = torch.log_softmax(
                            model(batch_tokens_masked.cuda())["logits"], dim=-1
                        )
                    all_token_probs.append(token_probs[:, i])  # vocab size
                token_probs = torch.cat(all_token_probs, dim=0).unsqueeze(0)
                df[model_location] = df.apply(
                    lambda row: label_row(
                        row[args.mutation_col],
                        seq,
                        token_probs,
                        alphabet,
                        offset,
                    ),
                    axis=1,
                )
            elif args.scoring_strategy == "pseudo-ppl":
                tqdm.pandas()
                df[model_location] = df.progress_apply(
                    lambda row: compute_pppl(
                        row[args.mutation_col], seq, model, alphabet, offset
                    ),
                    axis=1,
                )

    return df


def get_gene(gene, models):
    total_len = len(pep_dict[gene].seq)
    if total_len > 1000:
        # warnings.warn(f"{gene} is too long ({total_len} aa), spliting in to 500 aa chunks")
        dfs = []
        for i in range(total_len//500+1):
            offset = i*500+1
            seq =  str(pep_dict[gene].seq[i*500:min((i+1)*500,total_len)])
            df = get_saturated_mutations(seq, gene, offset)
            if len(seq) == 0:
                continue
            dfs.append(get_output(df, seq, offset, models, args))
        dfs = pd.concat(dfs).reset_index(drop=True)
        # dfs.to_feather(f"{args.out_dir}/{gene}.feather")
        return dfs
    else:
        seq = str(pep_dict[gene].seq)
        offset=1
        df = get_saturated_mutations(seq, gene, offset)
        dfs = get_output(df, seq, offset, models, args).reset_index(drop=True)
        # dfs.to_feather(f"{args.out_dir}/{gene}.feather")
        return dfs


def get_gene_with_fasta(gene, fasta, models):
    seqf = SeqIO.read(fasta, "fasta")
    total_len = len(seqf.seq)
    if total_len > 1000:
        warnings.warn(f"{gene} is too long ({total_len} aa), spliting in to 500 aa chunks")
        dfs = []
        for i in range(total_len//500+1):
            offset = i*500+1
            seq =  str(seqf.seq[i*500:min((i+1)*500,total_len)])
            df = get_saturated_mutations(seq, gene, offset)
            dfs.append(get_output(df, seq, offset, models, args))
        dfs = pd.concat(dfs).reset_index(drop=True)
        # dfs.to_feather(f"{args.out_dir}/{gene}.feather")
        return dfs
    else:
        seq = str(seqf.seq)
        offset=1
        df = get_saturated_mutations(seq, gene, offset)
        dfs = get_output(df, seq, offset, models, args).reset_index(drop=True)
        # dfs.to_feather(f"{args.out_dir}/{gene}.feather")
        return dfs


def get_saturated_mutations(seq, gene, offset=1):
    mut_list = []
    for i,a in enumerate(seq):
        for j in amino_acid_list:
            mut_list.append(a+str(i+offset)+j)
    data = pd.DataFrame({'Mutation':mut_list, 'Gene':gene})
    return data


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    uniprot_path = "/pmglocal/alb2281/es_paper/data/UP000005640_9606.fasta.gz"
    amino_acid_list = 'ACDEFGHIKLMNPQRSTVWY'
    model_name = args.model_location[0].split("/")[-1].split(".pt")[0]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if args.gene is not None:
        genes = [args.gene]
    else:
        with open(args.genes_file, "r") as f:
            data = f.read()
        genes = data.strip("\n").split("\n")
    with gzip.open(uniprot_path, 'rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_uniprot_accession
        )
    models = []

    for model_location in args.model_location:
        model, alphabet = pretrained.load_model_and_alphabet(model_location)
        model.eval()
        if torch.cuda.is_available() and not args.nogpu:
            model = model.cuda()
            print("Transferred model to GPU...")

        batch_converter = alphabet.get_batch_converter()
        models.append((model, batch_converter, alphabet, model_location))

    all_gene_dfs = []
    for i in tqdm(genes):
        if args.fasta:
            df = get_gene_with_fasta(i, args.fasta, models)
        else:
            df = get_gene(i, models)

        all_gene_dfs.append(df.iloc[:, 0:3])
    
    all_gene_dfs = pd.concat(all_gene_dfs)
    all_gene_dfs.reset_index(drop=True).to_feather(f"{output_dir}/{model_name}.feather")

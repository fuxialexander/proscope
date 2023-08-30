# %%
import gzip
import json
import os
from collections import defaultdict
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from tqdm import tqdm

from proscope.protein import Protein, setup_global_variables

# %%
setup_global_variables()

# %%
def parse_atm_record(line):
    """Get the atm record"""
    record = defaultdict()
    record["name"] = line[0:6].strip()
    record["atm_no"] = int(line[6:11])
    record["atm_name"] = line[12:16].strip()
    record["atm_alt"] = line[17]
    record["res_name"] = line[17:20].strip()
    record["chain"] = line[21]
    record["res_no"] = int(line[22:26])
    record["insert"] = line[26].strip()
    record["resid"] = line[22:29]
    record["x"] = float(line[30:38])
    record["y"] = float(line[38:46])
    record["z"] = float(line[46:54])
    record["occ"] = float(line[54:60])
    record["B"] = float(line[60:66])

    return record


def read_pdb(pdbfile):
    """Read a pdb file predicted with AF and rewritten to conatin all chains"""

    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, "r") as file:
        for line in file:
            if not line.startswith("ATOM"):
                continue
            record = parse_atm_record(line)
            # Get CB - CA for GLY
            if record["atm_name"] == "CB" or (
                record["atm_name"] == "CA" and record["res_name"] == "GLY"
            ):
                if record["chain"] in [*chain_coords.keys()]:
                    chain_coords[record["chain"]].append(
                        [record["x"], record["y"], record["z"]]
                    )
                    chain_plddt[record["chain"]].append(record["B"])
                else:
                    chain_coords[record["chain"]] = [
                        [record["x"], record["y"], record["z"]]
                    ]
                    chain_plddt[record["chain"]] = [record["B"]]

    # Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return chain_coords, chain_plddt


def calc_pdockq(chain_coords, chain_plddt, t):
    """Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    """

    # Get coords and plddt per chain
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    # Calc 2-norm
    mat = np.append(coords1, coords2, axis=0)
    a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
    dists = np.sqrt(np.sum(a_min_b.T**2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1, l1:]  # upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists <= t)

    if contacts.shape[0] < 1:
        pdockq = 0
        ppv = 0
    else:
        # Get the average interface plDDT
        avg_if_plddt = np.average(
            np.concatenate(
                [plddt1[np.unique(contacts[:, 0])], plddt2[np.unique(contacts[:, 1])]]
            )
        )
        # Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt * np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018

        # PPV
        PPV = np.array(
            [
                0.98128027,
                0.96322524,
                0.95333044,
                0.9400192,
                0.93172991,
                0.92420274,
                0.91629946,
                0.90952562,
                0.90043139,
                0.8919553,
                0.88570037,
                0.87822061,
                0.87116417,
                0.86040801,
                0.85453785,
                0.84294946,
                0.83367787,
                0.82238224,
                0.81190228,
                0.80223507,
                0.78549007,
                0.77766077,
                0.75941223,
                0.74006263,
                0.73044282,
                0.71391784,
                0.70615739,
                0.68635536,
                0.66728511,
                0.63555449,
                0.55890174,
            ]
        )

        pdockq_thresholds = np.array(
            [
                0.67333079,
                0.65666073,
                0.63254566,
                0.62604391,
                0.60150931,
                0.58313803,
                0.5647381,
                0.54122438,
                0.52314392,
                0.49659878,
                0.4774676,
                0.44661346,
                0.42628389,
                0.39990988,
                0.38479715,
                0.3649393,
                0.34526004,
                0.3262589,
                0.31475668,
                0.29750023,
                0.26673725,
                0.24561247,
                0.21882689,
                0.19651314,
                0.17606258,
                0.15398168,
                0.13927677,
                0.12024131,
                0.09996019,
                0.06968505,
                0.02946438,
            ]
        )
        inds = np.argwhere(pdockq_thresholds >= pdockq)
        if len(inds) > 0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]

    return pdockq, ppv


class AFScore(object):
    """AlphaFold scores"""

    def __init__(self, chain_len, js_file) -> None:
        with open(js_file) as f:
            score = json.load(f)
        js_file = os.path.basename(js_file)
        self.name = js_file.split("scores")[0]
        self.model_type = js_file.split("rank_")[1][4:].split("_model")[0]
        self.model_num = int(
            js_file.split("rank_")[1][4:].split("_model_")[1].split("_seed")[0]
        )
        self.seed = int(
            js_file.split("rank_")[1][4:]
            .split("_model_")[1]
            .split("_seed_")[1]
            .split(".json")[0]
        )
        self.rank = int(js_file.split("rank_")[1].split("_")[0])
        # self.score = score
        self.chain_len = chain_len
        self.max_pae = score["max_pae"]
        self.iptm = score["iptm"]
        self.plddt = np.array(score["plddt"])
        self.pae = np.array(score["pae"])
        self.ptm = np.array(score["ptm"])
        self.interchain_min_pae = self._get_interchain_min_pae()

    def _get_interchain_min_pae(self):
        # first set all non-interchain pae to 100, then get the min
        interchain_pae = self.pae.copy()
        # set first chain to 100
        interchain_pae[: self.chain_len[0], : self.chain_len[0]] = 100
        for i in range(len(self.chain_len)):
            start = sum(self.chain_len[:i])
            end = sum(self.chain_len[: i + 1])
            interchain_pae[start:end, start:end] = 100
        return np.min(interchain_pae)


# %%
class AFResult(object):
    """AlphaFold result from a files"""

    def __init__(self, result_dir, fasta_path) -> None:
        self.dir = result_dir
        self.fasta_path = fasta_path
        self.fasta = self._parse_fasta()

        self.name = self.fasta.id
        self.config = self._parse_config()
        self.pae = self._parse_pae()
        self.pdb = ""
        self.pdbs = os.path.join(
            self.dir,
            self.name
            + "_unrelaxed_rank_*_{model_type}_model_*.pdb".format(
                model_type=self.config["model_type"]
            ),
        )
        self.pdb = glob(self.pdbs)[0]
        pdockq_max = 0
        ppv_max = 0
        for f in glob(self.pdbs):
            chain_coords, chain_plddt = read_pdb(f)
            t = 8  # Distance threshold, set to 8 Å
            pdockq, ppv = calc_pdockq(chain_coords, chain_plddt, t)
            if pdockq > pdockq_max:
                pdockq_max = pdockq
                ppv_max = ppv
                self.pdb = f
        # get chain length from self.pdb
        chain_coords, chain_plddt = read_pdb(self.pdb)
        self.chain_len = [len(chain_coords[c]) for c in chain_coords]
        self.pdockq = pdockq_max
        self.ppv = ppv_max
        # self.fasta = self._parse_fasta()
        # self.chains = str(self.fasta.seq).split(':')
        # self.chain_len = [len(c) for c in self.chains]
        self.scores = self._parse_scores()
        self.interchain_min_pae = np.array(
            [s.interchain_min_pae for s in self.scores]
        ).min()
        self.max_pae = np.array([s.max_pae for s in self.scores]).max()
        self.min_pae = np.array([s.pae.min() for s in self.scores]).min()
        self.iptm = np.array([s.iptm for s in self.scores]).max()
        self.plddt = np.array([s.plddt for s in self.scores]).max(axis=0)
        self.mean_plddt = np.mean(self.plddt)

    def _parse_fasta(self):
        basename = os.path.basename(self.dir)
        if os.path.exists(self.fasta_path):
            fasta = SeqIO.read(self.fasta_path, "fasta")
        else:
            gene_seq = globals()['seq'][globals()['genename_to_uniprot'][basename]]
            fasta = SeqIO.SeqRecord(seq=gene_seq + ":" + gene_seq, id=basename)
        return fasta

    def _parse_config(self):
        with open(os.path.join(self.dir, "config.json")) as f:
            config = json.load(f)
        return config

    def _parse_pae(self):
        with open(
            os.path.join(self.dir, self.name + "_predicted_aligned_error_v1.json")
        ) as f:
            pae = json.load(f)
        return np.array(pae["predicted_aligned_error"])

    def _parse_scores(self):
        scores = []
        file_pattern = os.path.join(
            self.dir,
            self.name
            + "_scores_rank_*_{model_type}_model_*.json".format(
                model_type=self.config["model_type"]
            ),
        )
        for js_file in glob(file_pattern):
            scores.append(AFScore(self.chain_len, js_file))
        return scores

    def __repr__(self) -> str:
        # print chain_len, mean plddt, iptm, max_pae, interchain_min_pae
        return "chain_len\t{chain_len}\nmean_plddt\t{mean_plddt}\niptm\t{iptm}\nmax_pae\t{max_pae}\nmin_pae\t{min_pae}\ni_min_pae\t{interchain_min_pae}\n".format(
            chain_len=self.chain_len,
            mean_plddt=self.mean_plddt,
            iptm=self.iptm,
            max_pae=self.max_pae,
            min_pae=self.min_pae,
            interchain_min_pae=self.interchain_min_pae,
        )


# %%
class AFMonomer(AFResult):
    """An AFResult class for AF2 Monomer prediction"""

    def __init__(self, result_dir, name, fasta_dir=None) -> None:
        super().__init__(result_dir, name, fasta_dir)


class AFHomodimer(AFResult):
    """An AFResult class for AF2 Homodimer prediction"""


class AFMultimer(AFResult):
    """An AFResult class for AF2 Multimer prediction"""


# %%
class AFPairseg(object):
    """
    Pass the base result directory name Gene1_Gene2. In the directory there are many subfolders like Gene1_1_Gene2_2.
    Each subfolder contains the AF2 result for the pair of segments between the two genes. The subfolder name is the name for AFResult.
    The segment corresponds to the low_and_high_plddt_region_sequence in the Protein class.
    """

    def __init__(self, results_root_dir, fasta_root_dir) -> None:
        self.results_root_dir = results_root_dir
        self.fasta_root_dir = fasta_root_dir
        self.name = os.path.basename(results_root_dir)
        self.gene1 = self.name.split("_")[0]
        self.gene2 = self.name.split("_")[1]
        self.protein1 = Protein(self.gene1)
        self.protein2 = Protein(self.gene2)
        self.pairs_score, self.pairs_data = self._parse_pairs()

    def _init_score(self, len1, len2):
        score = {}
        score["mean_plddt"] = np.zeros((len1, len2))
        score["max_pae"] = np.zeros((len1, len2))
        score["min_pae"] = np.zeros((len1, len2))
        score["iptm"] = np.zeros((len1, len2))
        score["interchain_min_pae"] = np.zeros((len1, len2))
        score["pdockq"] = np.zeros((len1, len2))
        score["ppv"] = np.zeros((len1, len2))
        score["plddt"] = pd.DataFrame()
        return score

    def _parse_pairs(self):
        pairs = {}
        score = self._init_score(
            len(self.protein1.low_or_high_plddt_region),
            len(self.protein2.low_or_high_plddt_region),
        )

        for pair_dir in glob(os.path.join(self.results_root_dir, "*")):
            pair_name = os.path.basename(pair_dir)
            seg1 = int(pair_name.split("_")[1])
            seg2 = int(pair_name.split("_")[3])
            range1 = self.protein1.low_or_high_plddt_region[seg1]
            range2 = self.protein2.low_or_high_plddt_region[seg2]
            pair_fasta = os.path.join(self.fasta_root_dir, pair_name + ".fasta")
            res = AFResult(pair_dir, pair_fasta)
            pairs[pair_name] = res
            score["mean_plddt"][seg1, seg2] = res.mean_plddt
            score["max_pae"][seg1, seg2] = res.max_pae
            score["min_pae"][seg1, seg2] = res.min_pae
            score["iptm"][seg1, seg2] = res.iptm
            score["interchain_min_pae"][seg1, seg2] = res.interchain_min_pae
            score["pdockq"][seg1, seg2] = res.pdockq
            score["ppv"][seg1, seg2] = res.ppv
            plddt_df = pd.DataFrame({"plddt": res.plddt})
            plddt_df["seg1"] = f"{self.gene1}_{str(seg1)}"
            plddt_df["seg2"] = f"{self.gene2}_{str(seg2)}"
            plddt_df["gene1"] = self.gene1
            plddt_df["gene2"] = self.gene2
            plddt_df["gene1_res"] = np.concatenate(
                [
                    np.arange(range1[0], range1[1]),
                    np.repeat(np.nan, range2[1] - range2[0]),
                ]
            )
            plddt_df["gene2_res"] = np.concatenate(
                [
                    np.repeat(np.nan, range1[1] - range1[0]),
                    np.arange(range2[0], range2[1]),
                ]
            )
            score["plddt"] = pd.concat([score["plddt"],plddt_df])
        return score, pairs

    def plot_plddt_gene1(self):
        fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        sns.lineplot(
            x=range(self.protein1.length),
            y=self.protein1.plddt,
            ax=ax,
            color="black",
            label="monomer",
            alpha=0.5,
            linewidth=10,
        )
        sns.lineplot(
            x="gene1_res",
            y="plddt",
            data=self.pairs_score["plddt"].query("~gene1_res.isna()"),
            hue="seg2",
            palette="tab20",
            ax=ax,
            linewidth=3,
        )
        # legend right outside, make line in legend thicker
        legend = ax.legend(bbox_to_anchor=(1.05, 1, 0.1, 0.1), loc=2, borderaxespad=0.0, fontsize=16)
        for handle in legend.legendHandles:
            handle.set_linewidth(5.0)
        # set xlabel to gene1 name
        ax.set_xlabel(self.gene1)
        ax.set_ylabel("pLDDT")
        plt.tight_layout()
        return fig, ax

    def plot_plddt_gene2(self):
        fig, ax = plt.subplots(1, 1, figsize=(20, 5))
        sns.lineplot(
            x=range(self.protein2.length),
            y=self.protein2.plddt,
            ax=ax,
            color="black",
            label="monomer",
            alpha=0.5,
            linewidth=10,
        )
        sns.lineplot(
            x="gene2_res",
            y="plddt",
            data=self.pairs_score["plddt"].query("~gene2_res.isna()"),
            hue="seg1",
            palette="tab20",
            ax=ax,
            linewidth=3,
        )
        # legend right outside, make line in legend thicker
        legend = ax.legend(bbox_to_anchor=(1.05, 1, 0.1, 0.1), loc=2, borderaxespad=0.0, fontsize=16)
        for handle in legend.legendHandles:
            handle.set_linewidth(5.0)
        # set xlabel to gene2 name
        ax.set_xlabel(self.gene2)
        ax.set_ylabel("pLDDT")
        plt.tight_layout()
        return fig, ax

    def plot_score_heatmap(self):
        fig, axs = plt.subplots(2, 2, figsize=(12, 12))
        sns.heatmap(
            self.pairs_score["interchain_min_pae"],
            cmap="Blues",
            ax=axs[0, 0],
            vmin=0,
            vmax=10,
        )
        sns.heatmap(
            self.pairs_score["mean_plddt"], cmap="Blues", ax=axs[0, 1], vmin=0, vmax=100
        )
        sns.heatmap(
            self.pairs_score["iptm"], cmap="Blues", ax=axs[1, 0], vmin=0, vmax=1
        )
        sns.heatmap(
            self.pairs_score["pdockq"], cmap="Blues", ax=axs[1, 1], vmin=0, vmax=0.5
        )

        axs[0, 0].set_title("interchain_min_pae")
        axs[0, 1].set_title("mean_plddt")
        axs[1, 0].set_title("iptm")
        axs[1, 1].set_title("pdockq")
        # set y label to gene1 and x to gene2
        for ax in axs.flat:
            ax.set(xlabel=f"{self.gene2}", ylabel=f"{self.gene1}")
        plt.tight_layout()
        return fig, axs

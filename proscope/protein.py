import gzip
import os
import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import xmlschema
from Bio import SeqIO
from Bio.PDB import PDBParser
from genericpath import exists
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from plotly import graph_objects as go

pio.templates.default = "plotly_white"
import matplotlib.pyplot as plt
import numpy as np

# gaussian filter
from scipy.ndimage import gaussian_filter1d
from tqdm import tqdm

bioparser = PDBParser()
from proscope.data import get_genename_to_uniprot, get_lddt, get_schema, get_seq

seq = get_seq()
genename_to_uniprot = get_genename_to_uniprot()
lddt = get_lddt()
schema = get_schema()


def generate_pair_sequence(P1, P2, output_dir):
    """generate pair sequence from row"""
    protein_a = Protein(P1)
    protein_b = Protein(P2)
    low_or_high_plddt_region_sequence_a = protein_a.low_or_high_plddt_region_sequence
    low_or_high_plddt_region_sequence_b = protein_b.low_or_high_plddt_region_sequence
    for i, seq_a in enumerate(low_or_high_plddt_region_sequence_a):
        for j, seq_b in enumerate(low_or_high_plddt_region_sequence_b):
            os.makedirs(f"{output_dir}/", exist_ok=True)
            filename = f"{output_dir}/{protein_a.gene_name}_{i}_{protein_b.gene_name}_{j}.fasta"
            with open(filename, "w") as f:
                f.write(
                    f">{protein_a.gene_name}_{i}.{protein_b.gene_name}_{j}\n{str(seq_a.seq)}:{str(seq_b.seq)}\n"
                )


def min_max(arr):
    """normalize an array to 0-1"""
    if isinstance(arr, np.ndarray):
        return (arr - arr.min()) / (arr.max() - arr.min())
    elif isinstance(arr, list):
        return [(a - min(arr)) / (max(arr) - min(arr)) for a in arr]
    elif isinstance(arr, pd.DataFrame):
        arr_copy = arr.copy()
        arr_copy['plddt'] = (arr_copy['plddt'] - arr_copy['plddt'].min()) / (arr_copy['plddt'].max() - arr_copy['plddt'].min())
        return arr_copy


def normalize(x):
    new_x = x - x.min()
    return new_x / new_x.max()


def smooth(x, window_size=10):
    result = gaussian_filter1d(x, sigma=window_size / 2)
    return normalize(result)


def get_3d_avg(x, pairwise_interaction):
    x = x * pairwise_interaction
    x = x.sum(1) / ((x > 0).sum(1) + 0.01)
    return x / x.max()


def square_grad(f):
    return normalize(np.gradient(f) ** 2)


def extract_wt_from_mut(str):
    return str[0:1]


def extract_alt_from_mut(str):
    return str[-1:]


def extract_pos_from_mut(str):
    return int(str[1:-1])


class Protein(object):
    """Protein class"""

    def __init__(self, gene_name, homodimer=False, use_es=False, 
                 esm_folder = "/manitou/pmg/users/xf2217/demo_data/esm1b/esm1b_t33_650M_UR90S_1",
                 af2_folder = "/manitou/pmg/users/xf2217/demo_data/af2",
                 window_size=10):
        """
        Args:
            gene_name (str): gene name
            homodimer (bool): whether to use homodimer structure
            use_es (bool): whether to use ES score
        """
        self.gene_name = gene_name
        self.uniprot_id = genename_to_uniprot[gene_name]
        self.plddt = lddt[self.uniprot_id]
        self.length = len(self.plddt)
        self.sequence = seq[self.uniprot_id]
        self.smoothed_plddt_gaussian = self.get_smoothed_plddt_gaussian()
        self.smoothed_plddt = self.get_smooth_plddt()
        self.domains = self.get_domain_from_uniprot()
        self.window_size = window_size
        self.esm_folder = esm_folder
        self.af2_folder = af2_folder

        if use_es:
            self.grad = self.get_plddt_grad()
            self.pairwise_distance = self.get_pairwise_distance(dimer=False, af2_folder=self.af2_folder)
            self.esm = self.get_esm()
            self.es = smooth(
                self.get_final_score_gated_grad_3d(), window_size=self.window_size
            )

    def get_domain_from_uniprot(self):
        """Get domain information from uniprot"""
        uniprot_id = self.uniprot_id
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        entry_dict = schema.to_dict(url)
        features = entry_dict["entry"][0]["feature"]
        df = []
        for feature in features:
            feature_type = feature["@type"]
            if feature_type == "chain":
                continue
            if "begin" in feature["location"]:
                if "@description" in feature:
                    feature_description = feature["@description"]
                else:
                    feature_description = feature["@type"]
                feature_begin = feature["location"]["begin"]["@position"]
                feature_end = feature["location"]["end"]["@position"]
            else:
                continue

            df.append(
                {
                    "feature_type": feature_type,
                    "feature_description": feature_description,
                    "feature_begin": feature_begin,
                    "feature_end": feature_end,
                }
            )

        return pd.DataFrame(df)

    def get_smooth_plddt(self, window_size=10):
        result = np.convolve(
            self.plddt, np.ones(window_size) / window_size, mode="same"
        )
        return result / np.max(result)

    def get_smoothed_plddt_gaussian(self, sigma=2):
        result = gaussian_filter1d(self.plddt, sigma=sigma)
        return result / np.max(result)

    def get_plddt_grad(self):
        grad = square_grad(self.smoothed_plddt_gaussian)
        grad = np.clip(grad, np.quantile(grad, 0.2), np.quantile(grad, 0.8))
        return grad

    def get_esm(
        self,
        format="pos",
    ):
        """Use UNIPROT_ID to get ESM score from ESM1b model"""
        # read PAX5 (PAX5) | Q02548.csv
        # variant,score,pos
        # M1K,-8.983,1
        # M1R,-8.712,1
        df = pd.read_csv(
            f"{self.esm_folder}/{self.uniprot_id}_LLR.csv",
            index_col=0,
        )
        if format == "long":
            # melt to long format, column, row, value
            df = df.reset_index().melt(id_vars="index")
            df["variant"] = df["variable"].str.replace(" ", "") + df["index"].astype(
                str
            )
            df["pos"] = df["variable"].apply(lambda x: int(x.split(" ")[1]))
            df = df.rename({"value": "esm"}, axis=1)
            df = df[["variant", "pos", "esm"]]
            return df
        elif format == "wide":
            return df
        elif format == "pos":
            return normalize(-df.mean(0).values)
        # df['ALT'] = df.variant.apply(extract_alt_from_mut)
        # df['REF'] = df.variant.apply(extract_wt_from_mut)
        # df['esm'] = df['esm'].astype(float)
        # df['pos'] = df['pos'].astype(int)
        # df = df[['esm', 'ALT', 'pos']].pivot_table(index='ALT', columns='pos', values='esm').fillna(0)
        # df = normalize(-df.mean(0).values)
        # return df

    @property
    def low_plddt_region(self, threshold=0.6):
        idx = np.where(self.smoothed_plddt < threshold)[0]
        # get regions from idx, join two regions if they are close (<30aa apart)
        regions = []
        for i in idx:
            if len(regions) == 0:
                regions.append([i, i + 1])
            else:
                if i - regions[-1][1] < 30:
                    regions[-1][1] = i
                else:
                    regions.append([i, i + 1])

        for i in regions:
            if i[1] - i[0] < 30:
                regions.remove(i)
        return regions

    @property
    def low_plddt_region_sequence(self, threshold=0.6):
        regions = self.low_plddt_region
        sequences = []
        for i, region in enumerate(regions):
            s = SeqIO.SeqRecord(
                self.sequence[region[0] : region[1]],
                id=self.gene_name + "_" + str(i),
                name=self.gene_name,
                description=self.gene_name
                + "_"
                + str(i)
                + ": "
                + str(region[0])
                + "-"
                + str(region[1]),
            )
            sequences.append(s)
        return sequences

    @property
    def low_or_high_plddt_region(self, threshold=0.6):
        """Get regions with low plddt and high plddt, define breakpoint by merge regions if they are close (<30bp apart)
        Note: if the last region is close to the end, merge it with the second last region
        If the second last region has already been removed, ignore it"""
        # get regions list from low_plddt_region, keep also high plddt regions
        region_breakpoint = list(np.concatenate(self.low_plddt_region))
        # add 0 if not in the list
        if region_breakpoint[0] != 0:
            region_breakpoint = [0] + region_breakpoint
        if region_breakpoint[-1] != self.length:
            region_breakpoint.append(self.length)
        # merge regions if they are close (<30bp apart)
        # special case: if the last region is close to the end, merge it with the second last region
        # (0, 15, 60, 200, 240, 264)
        len_breakpoint = len(region_breakpoint)
        region_breakpoint_output = region_breakpoint.copy()
        for i in range(len_breakpoint - 1):
            # len(region_breakpoint) = 6; range(5) = 0, 1, 2, 3, 4
            if region_breakpoint[i + 1] - region_breakpoint[i] < 30:
                if i == len_breakpoint - 2:  # i = 4, second last region
                    if region_breakpoint[i] in region_breakpoint_output:
                        region_breakpoint_output.remove(region_breakpoint[i])
                else:
                    if region_breakpoint[i + 1] in region_breakpoint_output:
                        region_breakpoint_output.remove(region_breakpoint[i + 1])
        regions = []
        for i in range(len(region_breakpoint_output) - 1):
            regions.append(
                [region_breakpoint_output[i], region_breakpoint_output[i + 1]]
            )

        return regions

    @property
    def low_or_high_plddt_region_sequence(self, threshold=0.6):
        # get regions list from low_plddt_region, keep also high plddt regions
        sequences = []
        for i, region in enumerate(self.low_or_high_plddt_region):
            s = SeqIO.SeqRecord(
                self.sequence[region[0] : region[1]],
                id=self.gene_name + "_" + str(i),
                name=self.gene_name,
                description=self.gene_name
                + "_"
                + str(i)
                + ": "
                + str(region[0])
                + "-"
                + str(region[1]),
            )
            sequences.append(s)
        return sequences

    def get_final_score_gated_grad_3d(self, interaction_threshold=20):
        f = self.grad * self.esm
        pairwise_interaction = self.pairwise_distance < interaction_threshold
        f[(self.smoothed_plddt < 0.5)] = 0
        f = get_3d_avg(f, pairwise_interaction)
        f = normalize(f)

        return f

    def get_pairwise_distance(
        self, dimer=False, af2_folder="/manitou/pmg/users/xf2217/demo_data/af2"
    ):
        if dimer:
            structure = bioparser.get_structure(
                "homodimer", f"{af2_folder}/dimer_structures/" + self.gene_name + ".pdb"
            )
            model = structure[0]
            chain = model["A"]
            residues = [r for r in model.get_residues()]
            whole_len = len(residues)
            chain_len = len(chain)
            distance = np.zeros((whole_len, whole_len))
            for i, residue1 in enumerate(residues):
                for j, residue2 in enumerate(residues):
                    # compute distance between CA atoms
                    try:
                        d = residue1["CA"] - residue2["CA"]
                        distance[i][j] = d
                        distance[j][i] = d
                    except KeyError:
                        continue
            distance = np.fmin(
                distance[0:chain_len, 0:chain_len],
                distance[0:chain_len, chain_len:whole_len],
            )

        else:
            if exists(f"{af2_folder}/pairwise_interaction/" + self.uniprot_id + ".npy"):
                distance = np.load(
                    f"{af2_folder}/pairwise_interaction/" + self.uniprot_id + ".npy"
                )
            else:
                # make sure structures folder exists
                if not exists(f"{af2_folder}/structures/"):
                    os.makedirs(f"{af2_folder}/structures/")

                # download pdb file from AFDB to structures/
                import urllib.request

                url = (
                    "https://alphafold.ebi.ac.uk/files/AF-"
                    + self.uniprot_id
                    + "-F1-model_v4.pdb"
                )
                urllib.request.urlretrieve(
                    url,
                    f"{af2_folder}/structures/AF-"
                    + self.uniprot_id
                    + "-F1-model_v4.pdb",
                )

                # https://alphafold.ebi.ac.uk/files/AF-Q02548-F1-model_v4.pdb
                with open(
                    f"{af2_folder}/structures/AF-"
                    + self.uniprot_id
                    + "-F1-model_v4.pdb",
                    "r",
                ) as f:
                    structure = bioparser.get_structure("monomer", f)

                model = structure[0]
                chain = model["A"]
                residues = [r for r in model.get_residues()]
                whole_len = len(residues)
                chain_len = len(chain)
                distance = np.zeros((whole_len, whole_len))
                for i, residue1 in enumerate(residues):
                    for j, residue2 in enumerate(residues):
                        # compute distance between CA atoms
                        try:
                            d = residue1["CA"] - residue2["CA"]
                            distance[i][j] = d
                            distance[j][i] = d
                        except KeyError:
                            continue
            # make sure pairwise_interaction folder exists
            if not exists(f"{af2_folder}/pairwise_interaction/"):
                os.makedirs(f"{af2_folder}/pairwise_interaction/")

            np.save(
                f"{af2_folder}/pairwise_interaction/" + self.uniprot_id + ".npy",
                distance,
            )
        return distance

    def plotly_plddt(
        self,
        pos_to_highlight=None,
        to_compare=None,
        filename=None,
        show_low_plddt=True,
        show_domain=True,
        domains_to_show=[
            "region of interest",
            "DNA-binding region",
            "short sequence motif",
        ],
    ):
        # Initialize Plotly Figure
        fig = go.Figure()

        # Highlight low pLDDT regions
        if show_low_plddt:
            for i, region in enumerate(self.low_or_high_plddt_region):
                segment_name = f"{self.gene_name}_{str(i)} pLDDT={np.mean(self.smoothed_plddt[region[0]:region[1]]):.2f}"
                fig.add_trace(
                        go.Scatter(
                            x=[region[0], region[1], region[1], region[0]],
                            y=[0.8, 0.8, 1, 1],
                            fill='toself',
                            fillcolor='grey',
                            opacity=0.2,
                            line=dict(color='black'),
                            hoverinfo='text',
                            hovertext=[segment_name],
                            mode='lines',
                            showlegend=True,
                            legendgroup="pLDDT Segments",
                            legendgrouptitle=dict(text="pLDDT Segments"),
                            name=segment_name
                        )
                    )
        # Highlight specified positions
        if pos_to_highlight is not None:
            pos_to_highlight = np.array(pos_to_highlight) - 1
            fig.add_trace(
                go.Scatter(
                    x=pos_to_highlight,
                    y=self.smoothed_plddt[pos_to_highlight],
                    mode="markers",
                    marker=dict(color="orange", size=8),
                )
            )

        # Show domains if applicable, color by feature_type
        if show_domain:
            # Create a color mapping for each unique feature_type
            feature_types = self.domains.query('feature_type.isin(@domains_to_show)').feature_description.unique()
            # use tab20 color map
            colormap = plt.get_cmap("Set3").colors
            colors = colormap[1 : len(feature_types) + 1]
            # convert to rgb string
            colors = ["rgb" + str(i) for i in colors]
            color_mapping = dict(zip(feature_types, colors))

            for i, domain in self.domains.iterrows():
                if domain.feature_type in domains_to_show:
                    # Get the color for this feature_type
                    color = color_mapping[domain.feature_description]
                    # Add domain feature as filled rectangle
                    fig.add_trace(
                        go.Scatter(
                            x=[domain.feature_begin, domain.feature_end, domain.feature_end, domain.feature_begin],
                            y=[0.6, 0.6, 0.8, 0.8],
                            fill='toself',
                            fillcolor=color,
                            opacity=0.2,
                            line=dict(color=color),
                            hoverinfo='text',
                            hovertext=[domain.feature_description],
                            mode='lines',
                            showlegend=True,
                            legendgroup=domain.feature_type,
                            legendgrouptitle=dict(text=domain.feature_type),
                            name=domain.feature_description
                        )
                    )

        # Plot main pLDDT line
        fig.add_trace(
            go.Scatter(
                y=self.smoothed_plddt,
                mode="lines",
                name="pLDDT",
                line=dict(color="orange"),
            )
        )

        # Plot secondary comparison line if specified
        if to_compare is not None:
            if isinstance(to_compare, dict):
                for name, (color, data)  in to_compare.items():
                    fig.add_trace(
                        go.Scatter(
                            y=data,
                            mode="lines",
                            name=name,
                            line=dict(color=color),
                        )
                    )
            else:
                fig.add_trace(go.Scatter(y=to_compare, mode="lines", name="To Compare"))

        # Plot ES if available and no comparison data is specified
        elif hasattr(self, "es"):
            fig.add_trace(
                go.Scatter(y=self.es, mode="lines", name="ES", line=dict(color="blue"))
            )

        # Additional plot settings
        fig.update_layout(
            title=f"{self.gene_name} pLDDT",
            xaxis_title="Residue",
            yaxis_title="pLDDT",
            xaxis=dict(
                tickvals=np.arange(0, self.length, 100),
                ticktext=np.arange(1, self.length + 1, 100),
            ),
            font=dict(
                family="Arial",
            )
        )

        # Save figure if filename is provided
        if filename is not None:
            fig.write_image(filename)

        return fig

    def plot_plddt_manuscript(self, to_compare=None, filename=None):
        plt.figure(figsize=(12, 3))
        plt.plot(self.plddt)
        if to_compare is not None:
            plt.plot(to_compare)
        # highlight low plddt region
        for region in self.low_plddt_region:
            plt.axvspan(region[0], region[1], ymax=1, ymin=0, color="grey", alpha=0.1)

        # highlight domain, color by feature_type
        cmap = plt.get_cmap("tab20").colors
        # map feature_type to color
        feature_type_to_color = {}
        self.domains = self.domains.query(
            '(feature_type=="domain") or (feature_type=="region of interest")'
        )
        # for i, t in enumerate(obj.domains.feature_type.unique()):
        #     feature_type_to_color[t] = cmap[i+5]
        feature_type_to_color["domain"] = cmap[5]
        feature_type_to_color["region of interest"] = cmap[6]

        y_span = 0.1  # 0.8/len(obj.domains.feature_type.unique())
        for i, domain in self.domains.iterrows():
            idx = np.where(self.domains.feature_type.unique() == domain.feature_type)[
                0
            ][0]
            plt.axvspan(
                domain.feature_begin,
                domain.feature_end,
                ymax=idx * y_span + y_span,
                ymin=idx * y_span,
                color=feature_type_to_color[domain.feature_type],
                alpha=0.2,
            )
        # add legend of domain color
        legend_elements = []
        for i in self.domains.feature_type.unique():
            legend_elements.append(Patch(facecolor=feature_type_to_color[i], label=i))
            # reverse the order of legend
        legend_elements = legend_elements[::-1]
        # add "low plddt region" to legend
        legend_elements.append(Patch(facecolor="grey", label="low pLDDT region"))
        # add number index of low or high plddt region on top of the plot
        # for i, region in enumerate(obj.low_or_high_plddt_region):
        #     plt.text(region[0], 0.9, f"{i}", fontsize=12)
        # legend outside the plot
        plt.legend(
            handles=legend_elements,
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
        )
        plt.ylabel(f"{self.gene_name} pLDDT")
        plt.xlabel("Residue")
        # set xlim to the length of protein
        plt.xlim(0, len(self.plddt))
        # plt.ylabel("pLDDT")
        # plt.tight_layout()
        if filename is not None:
            plt.savefig(filename, dpi=300, bbox_inches="tight")
        plt.show()

    def plot_plddt(
        self,
        pos_to_highlight=None,
        to_compare=None,
        filename=None,
        show_low_plddt=True,
        show_domain=True,
        domains_to_show=["region of interest", "DNA-binding region", "splice variant"],
    ):
        fig, ax = plt.subplots(figsize=(20, 5))
        plt.plot(self.smoothed_plddt, label="pLDDT", color="orange")
        if to_compare is not None:
            plt.plot(to_compare)
        elif hasattr(self, "es"):
            plt.plot(self.es, label="ES", color="blue")
        if show_low_plddt:
            # highlight low plddt region
            for region in self.low_plddt_region:
                plt.axvspan(
                    region[0], region[1], ymax=1, ymin=0.8, color="red", alpha=0.2
                )

        if pos_to_highlight is not None:
            pos_to_highlight = np.array(pos_to_highlight) - 1
            if hasattr(self, "es"):
                plt.scatter(
                    pos_to_highlight, self.es[pos_to_highlight], color="blue", s=50
                )
            plt.scatter(
                pos_to_highlight,
                self.smoothed_plddt[pos_to_highlight],
                color="orange",
                s=50,
            )

        if show_domain:
            if domains_to_show is None:
                domains_to_show = self.domains.feature_type.unique()
            else:
                domains_to_show = np.array(domains_to_show)

            domains = self.domains.query("feature_type.isin(@domains_to_show)")
            # highlight domain, color by feature_type
            cmap = plt.get_cmap("Set3").colors
            # map feature_type to color
            feature_type_to_color = {}
            for i, t in enumerate(domains_to_show):
                feature_type_to_color[t] = cmap[i]

            y_span = 0.8 / len(domains_to_show)

            for i, domain in domains.query(
                "feature_type.isin(@domains_to_show)"
            ).iterrows():
                idx = np.where(domains_to_show == domain.feature_type)[0][0]
                plt.axvspan(
                    domain.feature_begin - 1,
                    domain.feature_end - 1,
                    ymax=idx * y_span + y_span,
                    ymin=idx * y_span,
                    color=feature_type_to_color[domain.feature_type],
                    alpha=0.2,
                )
            # add legend of domain color
            legend_elements = []
            for i in domains_to_show:
                legend_elements.append(
                    Patch(facecolor=feature_type_to_color[i], label=i)
                )
                # reverse the order of legend
            legend_elements = legend_elements[::-1]
            # legend outside the plot, append line plot legend in the end
            legend_elements.append(plt.Line2D([0], [0], color="blue", label="ES"))
            legend_elements.append(plt.Line2D([0], [0], color="orange", label="pLDDT"))
            plt.legend(
                handles=legend_elements,
                bbox_to_anchor=(1.05, 1, 0.1, 0.1),
                loc="upper left",
                borderaxespad=0.0,
                fontsize=16,
            )
        # add number index of low or high plddt region on top of the plot
        for i, region in enumerate(self.low_or_high_plddt_region):
            plt.text(region[0], 0.9, f"{i}", fontsize=12)

        plt.title(f"{self.gene_name} pLDDT")
        plt.xlabel("Residue")
        plt.ylabel("pLDDT")
        # set x ticks to start from 1
        plt.xticks(np.arange(0, self.length, 100), np.arange(1, self.length + 1, 100))
        plt.tight_layout()
        if filename is not None:
            plt.savefig(filename, dpi=300, bbox_inches="tight")
        return fig, ax

### Alejandro adds

def write_gene_list(feather_file, gene_list_file):
    gene_list = pd.read_feather(feather_file)["uniprot_id"].tolist()
    gene_list = list(set(gene_list))
    
    with open(gene_list_file, "w") as f:
        json.dump(gene_list, f)


def get_pairwise_distance(
    uniprot_id,
    dimer=False,
    af2_folder="/pmglocal/alb2281/es_paper/pairwise_dist",
    error_file=None,
):
    if dimer:
        structure = bioparser.get_structure(
            "homodimer", f"{af2_folder}/dimer_structures/" + uniprot_id + ".pdb"
        )
        model = structure[0]
        chain = model["A"]
        residues = [r for r in model.get_residues()]
        whole_len = len(residues)
        chain_len = len(chain)
        distance = np.zeros((whole_len, whole_len))
        for i, residue1 in enumerate(residues):
            for j, residue2 in enumerate(residues):
                # compute distance between CA atoms
                try:
                    d = residue1["CA"] - residue2["CA"]
                    distance[i][j] = d
                    distance[j][i] = d
                except KeyError:
                    continue
        distance = np.fmin(
            distance[0:chain_len, 0:chain_len],
            distance[0:chain_len, chain_len:whole_len],
        )
    else:
        if exists(f"{af2_folder}/pairwise_interaction/" + uniprot_id + ".npy"):
            distance = np.load(
                f"{af2_folder}/pairwise_interaction/" + uniprot_id + ".npy"
            )
        else:
            # make sure structures folder exists
            if not exists(f"{af2_folder}/structures/"):
                os.makedirs(f"{af2_folder}/structures/")

            # download pdb file from AFDB to structures/
            import urllib.request

            url = (
                "https://alphafold.ebi.ac.uk/files/AF-"
                + uniprot_id
                + "-F1-model_v4.pdb"
            )
            try:
                urllib.request.urlretrieve(
                    url,
                    f"{af2_folder}/structures/AF-"
                    + uniprot_id
                    + "-F1-model_v4.pdb",
                )
            except Exception as e:
                err_f.write(f"{uniprot_id},{e}\n")
                return

            # https://alphafold.ebi.ac.uk/files/AF-Q02548-F1-model_v4.pdb
            with open(
                f"{af2_folder}/structures/AF-"
                + uniprot_id
                + "-F1-model_v4.pdb",
                "r",
            ) as f:
                structure = bioparser.get_structure("monomer", f)

            model = structure[0]
            chain = model["A"]
            residues = [r for r in model.get_residues()]
            whole_len = len(residues)
            chain_len = len(chain)
            distance = np.zeros((whole_len, whole_len))
            for i, residue1 in enumerate(residues):
                for j, residue2 in enumerate(residues):
                    # compute distance between CA atoms
                    try:
                        d = residue1["CA"] - residue2["CA"]
                        distance[i][j] = d
                        distance[j][i] = d
                    except KeyError:
                        continue
        # make sure pairwise_interaction folder exists
        if not exists(f"{af2_folder}/pairwise_interaction/"):
            os.makedirs(f"{af2_folder}/pairwise_interaction/")

        np.save(
            f"{af2_folder}/pairwise_interaction/" + uniprot_id + ".npy",
            distance,
        )
        return


if __name__=="__main__":
    gene_list_file = "/manitou/pmg/users/alb2281/data/pairwise_dist/genes.json"
    error_file = "/pmglocal/alb2281/es_paper/pairwise_dist/error.csv"

    # write_gene_list(feather_file, gene_list_file)

    with open(gene_list_file, "r") as f: 
        gene_list = json.load(f)

    gene_list = [item for item in gene_list if item is not None]
    gene_list = sorted(gene_list)
    err_f = open(error_file, "w")
    for uniprot_id in tqdm(gene_list):
        get_pairwise_distance(uniprot_id, error_file=err_f)
    err_f.close()

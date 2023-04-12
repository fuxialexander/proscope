# base class of a protein, plddt-based domain detection and visualization
#%%
from Bio import SeqIO
import numpy as np
from matplotlib.patches import Patch
from tqdm import tqdm
import gzip
from matplotlib import pyplot as plt
import pandas as pd
import xmlschema
from af2 import AFResult
#%%
schema = xmlschema.XMLSchema('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd')
seq = {}
with gzip.open("data/uniprot_sprot.fasta.gz", 'rt') as f:
    for record in tqdm(SeqIO.parse(f, "fasta")):
        id = record.id.split('|')[1]
        seq[id] = record.seq

genename_to_uniprot = pd.read_csv(
        "data/uniprot_to_genename.txt", sep='\t').set_index('To').to_dict()['From']
        
lddt = dict()
with open("data/9606.pLDDT.tdt", 'r') as f:
    for line in f:
        id, score = line.strip().split('\t')
        lddt[id] = np.array(score.split(",")).astype(float)
#%%
class Protein(object):
    """Protein class"""
    def __init__(self, gene_name, homodimer=False):
        """
        Args:
        """
        self.gene_name = gene_name
        self.uniprot_id = genename_to_uniprot[gene_name]
        self.plddt = lddt[self.uniprot_id]
        self.length = len(self.plddt)
        # TODO: generalize the AFResult handling to all protein
        # if homodimer:
        #     dir = f"dimer_output/{gene_name}"
        #     result = AFResult(dir, gene_name)
        #     self.plddt = (result.plddt[0: self.length] + result.plddt[self.length: 2*self.length])/2
        self.sequence = seq[self.uniprot_id]
        self.smoothed_plddt = self.get_smooth_plddt()
        self.domains = self.get_domain_from_uniprot()
    
    def get_domain_from_uniprot(self):
        '''Get domain information from uniprot'''
        uniprot_id = self.uniprot_id
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        entry_dict = schema.to_dict(url)
        features = entry_dict['entry'][0]['feature']
        df = pd.DataFrame()
        for feature in features:
            feature_type = feature['@type']
            if feature_type == 'chain':
                continue
            if 'begin' in feature['location']:
                if '@description' in feature:
                    feature_description = feature['@description']
                else:
                    feature_description = feature['@type']
                feature_begin = feature['location']['begin']['@position']
                feature_end = feature['location']['end']['@position']
            else:
                continue

            df = df.append({
                'feature_type': feature_type,
                'feature_description': feature_description,
                'feature_begin': feature_begin,
                'feature_end': feature_end
            }, ignore_index=True)

        return df


    def get_smooth_plddt(self, window_size=10):
        result = np.convolve(self.plddt, np.ones(
            window_size)/window_size, mode='same')
        return result/np.max(result)


    def plot_plddt(self, to_compare=None):
        plt.figure(figsize=(20, 5))
        plt.plot(self.plddt)
        if to_compare is not None:
            plt.plot(to_compare)
        # highlight low plddt region
        for region in self.low_plddt_region:
            plt.axvspan(region[0], region[1], ymax=1, ymin=0.8, color='red', alpha=0.2)
        
        # highlight domain, color by feature_type
        cmap = plt.get_cmap('tab20').colors
        # map feature_type to color
        feature_type_to_color = {}
        for i, t in enumerate(self.domains.feature_type.unique()):
            feature_type_to_color[t] = cmap[i]

        y_span = 0.8/len(self.domains.feature_type.unique())
        for i, domain in self.domains.iterrows():
            idx = np.where(self.domains.feature_type.unique()==domain.feature_type)[0][0]
            plt.axvspan(domain.feature_begin, domain.feature_end, ymax=idx * y_span + y_span, ymin=idx * y_span 
             , color=feature_type_to_color[domain.feature_type], alpha=0.2)
        # add legend of domain color
        legend_elements = []
        for i in self.domains.feature_type.unique():
            legend_elements.append(Patch(facecolor=feature_type_to_color[i], label=i))
        # add number index of low or high plddt region on top of the plot
        for i, region in enumerate(self.low_or_high_plddt_region):
            plt.text(region[0], 0.9, f"{i}", fontsize=12)
        # legend outside the plot
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.title(f"{self.gene_name} pLDDT")
        plt.xlabel("Residue")
        plt.ylabel("pLDDT")

        plt.show()

    @property
    def low_plddt_region(self, threshold=0.6):
        idx = np.where(self.smoothed_plddt < threshold)[0]
        # get regions from idx, join two regions if they are close (<10bp apart) 
        regions = []
        for i in idx:
            if len(regions) == 0:
                regions.append([i, i+1])
            else:
                if i - regions[-1][1] < 30:
                    regions[-1][1] = i
                else:
                    regions.append([i, i+1])
        
        for i in regions:
            if i[1] - i[0] < 30:
                regions.remove(i)
        return regions

    @property
    def low_plddt_region_sequence(self, threshold=0.7):
        regions = self.low_plddt_region
        sequences = []
        for i, region in enumerate(regions):
            s = SeqIO.SeqRecord(self.sequence[region[0]:region[1]], id = self.gene_name + "_" + str(i), name=self.gene_name, description=self.gene_name+"_"+str(i)+": "+str(region[0])+"-"+str(region[1]))
            sequences.append(s)
        return sequences
    @property
    def low_or_high_plddt_region(self, threshold=0.7):
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
        for i in range(len_breakpoint-1):
            # len(region_breakpoint) = 6; range(5) = 0, 1, 2, 3, 4
            if region_breakpoint[i+1] - region_breakpoint[i] < 30:
                if i == len_breakpoint-2: # i = 4, second last region
                    if region_breakpoint[i] in region_breakpoint_output:
                        region_breakpoint_output.remove(region_breakpoint[i])
                else:
                    if region_breakpoint[i+1] in region_breakpoint_output:
                        region_breakpoint_output.remove(region_breakpoint[i+1])
        regions = []
        for i in range(len(region_breakpoint_output)-1):
            regions.append([region_breakpoint_output[i], region_breakpoint_output[i+1]])

        return regions

    @property
    def low_or_high_plddt_region_sequence(self, threshold=0.7):
        # get regions list from low_plddt_region, keep also high plddt regions
        sequences = []
        for i, region in enumerate(self.low_or_high_plddt_region):
            s = SeqIO.SeqRecord(self.sequence[region[0]:region[1]], id = self.gene_name + "_" + str(i), name=self.gene_name, description=self.gene_name+"_"+str(i)+": "+str(region[0])+"-"+str(region[1]))
            sequences.append(s)
        return sequences

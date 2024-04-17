from itertools import product
from ..sample import Sample
import numpy as np
import matplotlib.pyplot as plt
from .mutation_preprocessing import get_passed_SNVs
from pyfaidx import Fasta
import numpy as np
from collections import Counter, defaultdict
import logging
import tqdm
import subprocess
import os
import tempfile
logger = logging.getLogger('Length-DataReader')
logger.setLevel(logging.INFO)

LENGTH_BINS = np.array([0, 50, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 700])

LENGTH_INTERVALS = [f"[{LENGTH_BINS[i]}, {LENGTH_BINS[i+1]})]" for i in range(len(LENGTH_BINS)-1)]

nucleotide_order = ['A','C','G','T']

CONTEXTS = np.zeros(1)

CONTEXT_IDX = dict(zip(CONTEXTS, range(len(CONTEXTS))))

MUTATIONS = {
    context : LENGTH_INTERVALS
    for context in CONTEXTS
}

MUTATIONS_IDX = {
    context : dict(zip(m, range(len(MUTATIONS))))
    for context, m in MUTATIONS.items()
}






cmap = plt.colormaps['tab10']

_context_palette = {
    'A': cmap(0.3),
    'C': cmap(0.2),
    'G': cmap(0.9),
    'T': cmap(0.4)
}

ALPHA_LIST = [0.5 if fourmer[1] in ['A','G'] else 1 for fourmer in CONTEXTS]
COLOR_LIST = [_context_palette[fourmer[0]] for fourmer in CONTEXTS]
for i in range(len(COLOR_LIST)):
    color_tuple = list(COLOR_LIST[i])
    color_tuple[-1] = ALPHA_LIST[i]
    COLOR_LIST[i] = tuple(color_tuple)



class WeirdMutationError(Exception):
    pass


class LengthSample(Sample):

    N_CARDINALITY=1
    N_CONTEXTS=1
    N_MUTATIONS=len(LENGTH_BINS)
    N_ATTRIBUTES=1


    def plot(self, ax=None, figsize=(30,3), show_strand=True,**kwargs):
        #print(CONTEXTS)
        context_dist = np.zeros((self.N_CONTEXTS,))
        mutation_dist = np.zeros((self.N_CONTEXTS, self.N_MUTATIONS,))

        for context_idx, mutation_idx, weight in zip(
            self.context, self.mutation, self.weight
        ):
            context_dist[context_idx] += weight
            mutation_dist[context_idx, mutation_idx] += weight

        return self.plot_factorized(context_dist, mutation_dist, None, 
                                   ax=ax, 
                                   figsize=figsize,
                                   show_strand=show_strand, 
                                   **kwargs)
        

    @staticmethod
    def plot_factorized(context_dist, mutation_dist, attribute_dist, 
                        ax=None, figsize=(5,3), show_strand=True, fontsize=5, show_xticks=True, **kwargs):

        #joint_prob = (context_dist[:,None]*mutation_dist).ravel() # CxM
        #event_name = [(to_cosmic_str(c,m),'f') if c[1] in 'TC' else (to_cosmic_str(revcomp(c), complement[m]), 'r')
        #              for c in CONTEXTS for m in MUTATIONS[c]
        #             ]
        
        #event_prob = dict(zip(event_name, joint_prob))

        #fwd_events = np.array([event_prob[(event, 'f')] for event in COSMIC_SORT_ORDER])
        #rev_events = np.array([event_prob[(event, 'r')] for event in COSMIC_SORT_ORDER])

        
        if ax is None:
            _, ax = plt.subplots(1,1,figsize= figsize)

        plot_kw = dict(
            x = range(len(CONTEXTS)),
            color = COLOR_LIST,
            width = 1,
            edgecolor = 'white',
            linewidth = 0.5,
            #error_kw = {'alpha' : 0.5, 'linewidth' : 0.5}
        )
        extent = max(context_dist)

        
        ax.bar(height = context_dist, **plot_kw)
        ax.set(yticks = [0,extent], xticks = [], 
               xlim = (-1,len(CONTEXTS)), ylim = (-1e-6, 1.1*extent))
        if show_xticks:
            ax.set_xticks(range(len(CONTEXTS)), CONTEXTS, fontsize=fontsize)
        ax.tick_params(axis='x', labelrotation=90)
        ax.axhline(0, color = 'lightgrey', linewidth = 0.25)

        for s in ['left','right','top','bottom']:
            ax.spines[s].set_visible(False)

        return ax
    

    @classmethod
    def featurize_mutations(cls, 
                    motif_file, regions_file, fasta_file,
                    chr_prefix = '', 
                    weight_col = None, 
                    mutation_rate_file=None
                ):

        def process_line(line, fasta_object):
            
            fields = line.strip().split('\t')
            chrom=fields[0]; locus_idx=int(fields[3]); motif=fields[-1]
            frag_start=int(fields[-5]); frag_end=int(fields[-4]) 
            
            # Will be used later
            gc_content = float(fields[-2])

            length = frag_end - frag_start
            bin_index = np.digitize(length, LENGTH_BINS, right=False) - 1
            bin_index = min(LENGTH_BINS, len(LENGTH_INTERVALS)-1)
            return {
                'chrom' : chrom,
                'locus' : locus_idx,
                'mutation' : bin_index,
                'context' : 0,
                'attribute' : 0, # placeholder for now
                'weight' : 1, # currently no support for incorporating weights -- float(weight)
                'pos' : int(frag_start), # irrelevant since we merge mutations
                'cardinality': 0 # not using cardinality for lengths
            }
        


        query_statement = chr_prefix + '%CHROM\t%POS0\t%POS0\t%POS0|%REF|%ALT|' \
                                                    + ('1\n' if weight_col is None else f'%INFO/{weight_col}\n')

        try:
            # Keep the fixes from motif corpus ingestion and only select fragments that are fully contained in a particular locus
            
            # fix intersecting codes for fragment input (start position of fragment should be exist in each locus being intersected)
            # Column number in awk is offset by the number of columns in regions_file
            # Calculate the number of columns in regions_file
            with open(regions_file, 'r') as f:
                # Assuming that the file is not empty and has a header or at least one line
                num_cols = len(f.readline().strip().split())
            
            # Construct the awk command with the correct column number for the start/end position from fragment input file
            awk_cmd = f"awk '{{if (($2 <= ${num_cols + 2} && $2 <= ${num_cols + 3})&& (${num_cols + 2} < $3 && ${num_cols + 3} < $3)) print}}'" # both start and end pos of fragment should be exist in locus
            
            # Now construct the full command, incorporating the dynamically constructed awk_cmd
            cmd = (
                f"bedtools intersect -a {regions_file} -b - -sorted -wa -wb -split | "
                f"{awk_cmd}"
            )

            intersect_process = subprocess.Popen(
                cmd,
                shell=True, 
                stdin=open(motif_file),
                stdout=subprocess.PIPE,
                universal_newlines=True,
                bufsize=10000,
            )
            # intersect_process = subprocess.Popen(
            #     ['bedtools',
            #     'intersect',
            #     '-a', regions_file, 
            #     '-b', '-', 
            #     '-sorted',
            #     '-wa','-wb',
            #     '-split'],
            #     stdin=open(motif_file),
            #     stdout=subprocess.PIPE,
            #     universal_newlines=True,
            #     bufsize=10000,
            # )
            mutations_grouped = {}
            last_line =  ""
            max_locus_processed = 0
            with Fasta(fasta_file) as fa:

                while True:
                    line = intersect_process.stdout.readline()
                    if not line:
                        break
            
                    line_dict = process_line(line, fa)
                    last_line = line
                    max_locus_processed = max(max_locus_processed, int(line_dict['locus']))
                    mutation_group_key = f"{line_dict['chrom']}:{line_dict['locus']}:{line_dict['mutation']}"

                    if mutation_group_key not in mutations_grouped:
                        mutations_grouped[mutation_group_key] = {}

                    for key in line_dict:
                        if key == "weight":
                            mutations_grouped[mutation_group_key]["weight"] = mutations_grouped[mutation_group_key].get("weight", 0) + \
                                                                              line_dict["weight"]
                        else:
                            mutations_grouped[mutation_group_key][key] = line_dict[key] # so pos & attribute are constantly overwritten, can't be used as of now

            intersect_process.communicate()

            mutations = defaultdict(list)

            for key in mutations_grouped:
                mutation_dict = mutations_grouped[key]
                for k, v in mutation_dict.items():
                    mutations[k].append(v)
            
            for k, v in mutations.items():
                mutations[k] = np.array(v).astype(LengthSample.type_map[k])

            return cls(
                **mutations,
                name = os.path.abspath(motif_file),
                )
        
        finally:
            if not mutation_rate_file is None:
                clustered_vcf.close()


    @classmethod
    def get_context_frequencies(cls, window_set, fasta_file, n_jobs = 1):
        fournuc_matrix = np.zeros((len(window_set), 1, 1))
        # L, D, C
        fournuc_matrix = fournuc_matrix.transpose(((1,2,0))) # DON'T (!) add a pseudocount

        return fournuc_matrix # DON'T (!) normalize, the number of contexts in a window is part of the likelihood
    

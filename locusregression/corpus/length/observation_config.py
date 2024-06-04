from itertools import product
from ..sample import Sample
from ..reader_utils import read_windows
import numpy as np
import matplotlib.pyplot as plt
from pyfaidx import Fasta
from collections import Counter, defaultdict
import logging
import tqdm
import subprocess
import os
import tempfile

logger = logging.getLogger('Length-DataReader')
logger.setLevel(logging.INFO)

# Length bins
LENGTH_BINS = [(50, 100), (100, 105), (105, 110), (110, 115), (115, 120), 
               (120, 125), (125, 130), (130, 135), (135, 140), (140, 145), 
               (145, 150), (150, 155), (155, 160), (160, 165), (165, 170), 
               (170, 175), (175, 180), (180, 185), (185, 190), (190, 195), 
               (195, 200), (200, 210), (210, 220), (220, 230), (230, 240), 
               (240, 250), (250, 260), (260, 270), (270, 280), (280, 290), 
               (290, 300), (300, 310), (310, 320), (320, 330), (330, 340), 
               (340, 350), (350, 700)]

N_CARDINALITY = 1
N_CONTEXTS = 1
N_MUTATIONS = len(LENGTH_BINS)
N_ATTRIBUTES = 1

# Create mutation indices
MUTATION_IDX = {i: idx for idx, (i, j) in enumerate(LENGTH_BINS)}

class WeirdMutationError(Exception):
    pass

class LengthSample(Sample):
    N_CARDINALITY = N_CARDINALITY
    N_CONTEXTS = N_CONTEXTS
    N_MUTATIONS = N_MUTATIONS
    N_ATTRIBUTES = N_ATTRIBUTES

    def plot(self, ax=None, figsize=(30, 3), show_strand=True, **kwargs):
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
                        ax=None, figsize=(5, 3), show_strand=True, fontsize=5, show_xticks=True, **kwargs):
        
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            fig.set_dpi(300)

        plot_kw = dict(
            x=range(len(LENGTH_BINS)),
            color='tab:blue',
            width=1,
            edgecolor='white',
            linewidth=0.5,
        )
        extent = max(context_dist)

        ax.bar(height=context_dist, **plot_kw)
        ax.set(yticks=[0, extent], xticks=[], 
               xlim=(-1, len(LENGTH_BINS)), ylim=(-1e-6, 1.1 * extent))
        if show_xticks:
            ax.set_xticks(range(len(LENGTH_BINS)))
            ax.set_xticklabels([f'{i}-{j}' for i, j in LENGTH_BINS], fontsize=fontsize, rotation=90)
        else:
            ax.set_xticklabels([])
        ax.axhline(0, color='lightgrey', linewidth=0.25)

        for s in ['left', 'right', 'top', 'bottom']:
            ax.spines[s].set_visible(False)

        return ax
    

    @classmethod
    def featurize_mutations(cls, 
                            motif_file, regions_file, fasta_file,
                            chr_prefix='', 
                            weight_col=None, 
                            mutation_rate_file=None,
                            sample_weight=1.,
                            in_corpus=True,
                            **kw):
        
        def process_line(line, fasta_object, positive_file=True, in_corpus=True):
            fields = line.strip().split('\t')
            chrom = fields[0]
            locus_idx = int(fields[3])
            frag_start = int(fields[5])
            frag_end = int(fields[6])
            length = frag_end - frag_start

            for min_length, max_length in LENGTH_BINS:
                if min_length <= length < max_length:
                    mutation_idx = MUTATION_IDX[(min_length, max_length)]
                    break
            else:
                raise WeirdMutationError(f'Fragment length {length} not in defined bins.')

            return {
                'chrom': chrom,
                'locus': locus_idx,
                'mutation': mutation_idx,
                'context': 0,  # Only one context
                'attribute': 0,
                'weight': 1,
                'pos': frag_start,
                'cardinality': 0
            }

        num_cols = 4
        if in_corpus:
            awk_cmd = f"awk '{{if ($2 <= ${num_cols + 2} && $3 > ${num_cols + 2}) print}}'"
        else:
            awk_cmd = f"awk '{{if ($2 <= ${num_cols + 2}-1 && $3 > ${num_cols + 2}-1) print}}'"

        segments = []
        for region in read_windows(regions_file):
            for chr, start, end in region.segments():
                segments.append((chr, start, end, region.name))

        segments = sorted(segments, key=lambda x: (x[0], x[1]))

        with tempfile.NamedTemporaryFile() as temp_file:
            with open(temp_file.name, 'w') as output:
                for row in segments:
                    print(*row, sep='\t', file=output)

            cmd = (
                f"bedtools intersect -a {temp_file.name} -b - -sorted -wa -wb | "
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

            positive_file = ("pos_" in motif_file)
            mutations_grouped = {}
            max_locus_processed = 0
            with Fasta(fasta_file) as fa:
                while True:
                    line = intersect_process.stdout.readline()
                    if not line:
                        break
                    
                    line_dict = process_line(line, fa, positive_file=positive_file, in_corpus=in_corpus)
                    max_locus_processed = max(max_locus_processed, int(line_dict['locus']))
                    mutation_group_key = f"{line_dict['chrom']}:{line_dict['locus']}:{line_dict['context']}"

                    if mutation_group_key not in mutations_grouped:
                        mutations_grouped[mutation_group_key] = {}

                    for key in line_dict:
                        if key == "weight":
                            mutations_grouped[mutation_group_key]["weight"] = mutations_grouped[mutation_group_key].get("weight", 0) + line_dict["weight"]
                        else:
                            mutations_grouped[mutation_group_key][key] = line_dict[key]

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
            name=os.path.abspath(motif_file),
        )

    @classmethod
    def get_context_frequencies(cls, window_set, fasta_file, n_jobs=1, in_corpus=True):
        def count_contexts(bed12_region, fasta_object):
            # We always return 1 for the single context
            return [[1]]

        with Fasta(fasta_file) as fasta_object:
            context_matrix = [
                count_contexts(w, fasta_object) 
                for w in tqdm.tqdm(window_set, nrows=100, desc='Counting contexts')
            ]

        context_matrix = np.array(context_matrix).transpose((1, 2, 0))

        return context_matrix

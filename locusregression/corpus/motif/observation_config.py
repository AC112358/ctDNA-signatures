from itertools import product
from ..sample import Sample
from ..reader_utils import read_windows
import numpy as np
import matplotlib.pyplot as plt
from pyfaidx import Fasta
import numpy as np
from collections import Counter, defaultdict
import logging
import tqdm
import subprocess
import os
import tempfile
import pandas as pd
logger = logging.getLogger('Motif-DataReader')
logger.setLevel(logging.INFO)

trans_table = str.maketrans('ATCG', 'TAGC')
def comp(seq):
    return seq.translate(trans_table)
def revcomp(seq):
    return seq[::-1].translate(trans_table)


nucleotide_order = ['A','C','G','T']

CONTEXTS = sorted(
    map(lambda x : ''.join(x), product('ATCG','ATCG','ATCG', 'ATCG')), 
    key = lambda x : (x[0], x[1], x[2], x[3])
    )

CONTEXT_IDX = dict(zip(CONTEXTS, range(len(CONTEXTS))))

MUTATIONS = {
    context : [alt for alt in 'ACGT' if not alt == context[1]] # nucleotide cannot mutate to itself
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


class MotifSampleBase(Sample):

    N_CARDINALITY=1
    N_CONTEXTS=256
    N_MUTATIONS=1
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
            fig, ax = plt.subplots(1,1,figsize= figsize)
            fig.set_dpi(300)

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
            ax.set_xticks(range(len(CONTEXTS)))  # Set the ticks locations
            ax.set_xticklabels(CONTEXTS, fontsize=fontsize, rotation=90)  # Set the labels with rotation
        else:
            ax.set_xticklabels([])
        ax.axhline(0, color = 'lightgrey', linewidth = 0.25)

        for s in ['left','right','top','bottom']:
            ax.spines[s].set_visible(False)

        return ax
    

    @classmethod
    def featurize_mutations(cls, 
                    motif_file, regions_file, fasta_file,
                    **kw,
                ):

        def process_line(line, fasta_object, 
                         positive_file=True):
            
            fields = line.strip().split('\t')

            # @Sandra: in the future, these hard-coded indices will need to be replaced with a more robust method
            # Potentially, we could start from a BAM file instead of the custom fragment file format.
            # The first four columns will be chr,start,end,locus_idx, then the rest of the columns will be from the motif file.
            chrom=fields[0]; locus_idx=int(fields[3]); #motif=fields[-1]
            frag_start=int(fields[5]); frag_end=int(fields[6]) 
            
            '''# Will be used later
            #gc_content = float(fields[-2])

            in_4mer = ""
            out_4mer = ""

            if positive_file:
                out_4mer = motif[:4][::-1]
                in_4mer = motif[4:]
            else:
                out_4mer = comp(motif[4:])
                in_4mer = revcomp(motif[:4])
            
            context = ""
            if in_corpus:
                context = in_4mer
            else:
                context = out_4mer'''
            
            # @sandra: Can get the fragment end motifs from the fasta file - requires no preprocessing of fragment file
            if positive_file and cls.in_corpus:
                context = fasta_object[chrom][frag_start-1:frag_start+3].seq.upper()
            elif positive_file:
                context = fasta_object[chrom][frag_start-5 : frag_start-1][::-1].seq.upper()
            else:
                raise NotImplementedError("Only positive_file is supported")
            
            return {
                'chrom' : chrom,
                'locus' : locus_idx,
                'mutation' : 0,
                'context' : CONTEXT_IDX[context],
                'attribute' : 0, # placeholder for now
                'weight' : 1, # currently no support for incorporating weights -- float(weight)
                'pos' : int(frag_start), # irrelevant since we merge mutations
                'cardinality': 0 # not using cardinality for motifs
            }


        # fix intersecting codes for fragment input (start position of fragment should be exist in each locus being intersected)
        # Column number in awk is offset by the number of columns in regions_file
        # Calculate the number of columns in regions_file
        #with open(regions_file, 'r') as f:
        #    # Assuming that the file is not empty and has a header or at least one line
        #    num_cols = len(f.readline().strip().split())
        num_cols = 4
        # Construct the awk command with the correct column number for the start/end position from fragment input file

        if cls.in_corpus:
            # condition 1: start of region should match or be less than the start of the fragment
            # condition 2: end of region should be greater than the start of the fragment
            awk_cmd = f"awk '{{if ($2 <= ${num_cols + 2} && $3 > ${num_cols + 2}) print}}'" # for pos in5p, start pos of fragment should be exist in locus
        else:
            awk_cmd = f"awk '{{if ($2 <= ${num_cols + 2}-1 && $3 > ${num_cols + 2}-1) print}}'" # for pos out5p, start-1 pos of fragment should be exist in locus


        # @Sandra: changes to make mesoscale features work:
        # 1. Read the regions_file and break it into the "exon" segments
        # 2. Write each segment to a temporary bed file, but with the region name as the 4th column
        # 3. Use bedtools intersect to intersect the motif file with the temporary bed file
        # This allows you to run your awk command to check intersection rules, but for each exon segment instead of the whole region
        
        ##
        # 1.
        ## 
        segments = []
        for region in read_windows(regions_file):
            for chr, start, end in region.segments():
                segments.append((chr, start, end, region.name))

        segments = sorted(segments, key=lambda x: (x[0], x[1]))

        ##
        # 2.
        ##
        with tempfile.NamedTemporaryFile() as temp_file:
            with open(temp_file.name, 'w') as output:
                for row in segments:
                    print(*row, sep='\t', file=output)

            ##
            # 3.
            ##
            # Now construct the full command, incorporating the dynamically constructed awk_cmd
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
            # positive_file = True
            mutations_grouped = {}
            max_locus_processed = 0
            with Fasta(fasta_file) as fa:

                while True:
                    line = intersect_process.stdout.readline()
                    if not line:
                        break
            
                    line_dict = process_line(line, fa, positive_file=positive_file)
                    #last_line = line
                    max_locus_processed = max(max_locus_processed, int(line_dict['locus']))
                    mutation_group_key = f"{line_dict['chrom']}:{line_dict['locus']}:{line_dict['context']}"

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
            mutations[k] = np.array(v).astype(cls.type_map[k])

        return cls(
            **mutations,
            name = os.path.abspath(motif_file),
        )


    @classmethod
    def get_context_frequencies(cls, window_set, fasta_file, n_jobs = 1):
    
        def count_fournucs(bed12_region, fasta_object):

            def rolling(seq, w = 4):
                for i in range(len(seq) - w + 1):
                    yield seq[ i : i+w]

            fournuc_counts = Counter()
            N_counts=0
            # sandra fixed: for fragment input we don't need to consider only exonic region. so DON'T use bed12_region.segments()
            #chrom = bed12_region.chromosome
            #start = bed12_region.start
            #end = bed12_region.end
            
            for chrom, start, end in bed12_region.segments():
                if cls.in_corpus:
                    window_sequence = fasta_object[chrom][max(0, start-1) : end+3].seq.upper() # sandra fixed: for 4mer motifs it should be "[max(start-1,0) : end+3]" to avoid not counting the 4mer motifs starting after the end-3 index of the sequence. "start-1" is because of python indexing for list
                else:
                     window_sequence = fasta_object[chrom][max(start-1-3,0) : end+3].seq.upper()  # sandra: it becomes -3/+3 since we forced to include start-1 and end+1 in the locus for neg files

                for fournuc in rolling(window_sequence):
                    if not 'N' in fournuc:
                        fournuc_counts[fournuc]+=1

                    else:
                        N_counts+=1

                pseudocount = N_counts/(4**4)

            return [
                [fournuc_counts[context]+pseudocount for context in CONTEXTS]
            ] # sandra fixed: frequency is calculated in the 5' to 3' direction (not averaged with 3' to 5')
        
        with Fasta(fasta_file) as fasta_object:
            fournuc_matrix = [
                count_fournucs(w, fasta_object) 
                for w in tqdm.tqdm(window_set, nrows=100, desc = 'Aggregating 4mer content')
            ]

        # L, D, C
        fournuc_matrix = np.array(fournuc_matrix).transpose(((1,2,0))) # DON'T (!) add a pseudocount

        return fournuc_matrix # DON'T (!) normalize, the number of contexts in a window is part of the likelihood
    

class MotifSampleIn5p(MotifSampleBase):
    in_corpus = True
    @classmethod
    def get_context_frequencies(cls, window_set, fasta_file, n_jobs = 1):
        # for positive strand in5p motifs, no need to transform the sequence, but future negative strand motifs will need to be transformed
        return super(MotifSampleIn5p, cls).get_context_frequencies(window_set, fasta_file, n_jobs = n_jobs)


class MotifSampleOut5p(MotifSampleBase):
    in_corpus = False
    @classmethod
    def get_context_frequencies(cls, window_set, fasta_file, n_jobs = 1):
        transform_func = lambda x: x[::-1] # for positive strand out5p motifs
        _context_frequencies = super(MotifSampleOut5p, cls).get_context_frequencies(window_set, fasta_file,  n_jobs = n_jobs)
        context_frequencies = pd.DataFrame(_context_frequencies[0], index=CONTEXTS)
        context_frequencies.index = [transform_func(seq) for seq in context_frequencies.index]
        context_frequencies = context_frequencies.loc[CONTEXTS].values.reshape(_context_frequencies.shape)          
        return context_frequencies

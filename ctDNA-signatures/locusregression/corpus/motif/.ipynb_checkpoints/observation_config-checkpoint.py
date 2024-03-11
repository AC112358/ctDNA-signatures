from itertools import product
from ..sample import Sample
import numpy as np
import matplotlib.pyplot as plt

complement = {'A' : 'T','T' : 'A','G' : 'C','C' : 'G'}

def revcomp(seq):
    return ''.join(reversed([complement[nuc] for nuc in seq]))


def convert_to_cosmic(context, alt):

    if not context[1] in 'CT': 
        context, alt = revcomp(context), complement[alt]

    return context, alt


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

class MotifSample(Sample):

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
                        ax=None, figsize=(5,3), show_strand=True,**kwargs):

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
               xlim = (-1,len(CONTEXTS)), ylim = (-0.1, 1.1*extent))
        ax.set_xticks(range(len(CONTEXTS)), CONTEXTS)
        ax.tick_params(axis='x', labelrotation=90)
        ax.axhline(0, color = 'lightgrey', linewidth = 0.25)

        for s in ['left','right','top','bottom']:
            ax.spines[s].set_visible(False)

        return ax
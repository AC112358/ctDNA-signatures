import numpy as np
from abc import ABC, abstractmethod
from collections import defaultdict
import h5py as h5
import logging
from .sample import SampleLoader, InMemorySamples
from .sbs.observation_config import SBSSample
from .motif.observation_config import MotifSample
from pandas import DataFrame
from tqdm import trange
logger = logging.getLogger('Corpus')


class CorpusMixin(ABC):

    def __init__(self,
        metadata = {},
        exposures=None,*,
        type,
        name,
        samples,
        features,
        context_frequencies,
        regions,
    ):
        self.type = type
        self.name = name
        self.samples = samples
        self.features = features
        self.context_frequencies = context_frequencies
        self._shared_exposures = True
        self.metadata = metadata
        self.regions=regions
        
        if exposures is None:
            self._exposures = np.ones((1, self.locus_dim))
        else:
            assert exposures.shape == (1, self.locus_dim,)
            self._exposures = exposures

        assert context_frequencies.shape == \
            (self.cardinalities_dim, self.context_dim, self.locus_dim)
        

    @abstractmethod
    def observation_class(self):
        raise NotImplementedError()

    @property
    def shape(self):
        return {
            'context_dim' : self.context_dim,
            'mutation_dim' : self.mutation_dim,
            'attribute_dim' : self.attribute_dim,
            'locus_dim' : self.locus_dim,
            'feature_dim' : self.feature_dim,
            'cardinality_features_dim' : self.cardinality_features_dim,
            'cardinalities_dim' : self.cardinalities_dim,
        }
    
    @property
    def feature_names(self):
        return list(self.features.keys())

    @property
    def exposures(self):
        return self._exposures

    @property
    def shared_correlates(self):
        return True
    
    @property
    def shared_exposures(self):
        return self._shared_exposures
        
    @abstractmethod
    def __iter__(self):
        raise NotImplementedError()

    @abstractmethod
    def __len__(self):
        raise NotImplementedError()

    @abstractmethod
    def __getitem__(self, idx):
        raise NotImplementedError()

    @abstractmethod
    def subset_samples(self, idx):
        raise NotADirectoryError()

    @abstractmethod
    def subset_loci(self, loci):
        raise NotImplementedError()
    
    @abstractmethod
    def corpuses(self):
        raise NotImplementedError()
    
    @abstractmethod
    def corpus_names(self):
        raise NotImplementedError()
    
    @abstractmethod
    def get_corpus(self, name):
        raise NotImplementedError()


class Corpus(CorpusMixin):

    @classmethod
    def _get_observation_class(cls, classname):
        if classname.lower() == 'sbs':
            return SBSSample
        elif classname.lower() == 'fragment-motif':
            return MotifSample
        else:
            raise ValueError(f'Unknown corpus type {classname}')

    @property
    def observation_class(self):
        return self._get_observation_class(self.type)        
        
    @property
    def context_dim(self):
        return self.observation_class.N_CONTEXTS
    
    @property
    def mutation_dim(self):
        return self.observation_class.N_MUTATIONS
    
    @property
    def attribute_dim(self):
        return self.observation_class.N_ATTRIBUTES
    
    @property
    def locus_dim(self):
        return len(next(iter(self.features.values()))['values'])
    
    @property
    def feature_dim(self):
        return len([f for f in self.features.values() if not f['type'] == 'cardinality'])
    
    @property
    def cardinality_features_dim(self):
        return len([f for f in self.features.values() if f['type'] == 'cardinality'])
    
    @property
    def cardinalities_dim(self):
        return self.observation_class.N_CARDINALITY

    def __len__(self):
        return len(self.samples)


    def __getitem__(self, idx):
        sample = self.samples[idx]
        sample.corpus_name = self.name

        return sample


    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


    @property
    def corpuses(self):
        return [self]
    
    @property
    def corpus_names(self):
        return [self.name]
    
    @property
    def num_mutations(self):
        return sum([sum(sample.weight) for sample in self.samples])


    def subset_samples(self, subset_idx):

        return Corpus(
            type = self.type,
            samples = self.samples.subset(subset_idx),
            features=self.features,
            context_frequencies = self.context_frequencies,
            exposures = self.exposures,
            name = self.name,
            regions=self.regions,
        )


    def subset_loci(self, loci):

        subsample_lookup = dict(zip(loci, np.arange(len(loci)).astype(int)))
        
        bool_array = np.zeros(self.shape['locus_dim']).astype(bool)
        bool_array[loci] = True

        #total_mutations = 0
        new_samples = []
        for sample in self.samples:
                        
            mask = bool_array[sample.locus]

            new_sample = self.observation_class(**{
                'attribute' : sample.attribute[mask],
                'mutation' : sample.mutation[mask],
                'context' : sample.context[mask],
                'weight' : sample.weight[mask],
                'cardinality' : sample.cardinality[mask],
                'locus' : np.array([subsample_lookup[locus] for locus in sample.locus[mask]]).astype(int),
                'chrom' : sample.chrom[mask],
                'pos' : sample.pos[mask], 
                'name' : sample.name,
            })
            new_samples.append(new_sample)
        
        return Corpus(
            type = self.type,
            samples = InMemorySamples(new_samples),
            features = {
                feature_name : {
                    'type' : v['type'], 
                    'group' : v['group'],
                    'values' : v['values'][loci]
                }
                for feature_name, v in self.features.items()
            },
            context_frequencies = self.context_frequencies[:,:,loci],
            name = self.name,
            exposures = self.exposures[:, loci],
            regions=[self.regions[l] for l in loci]
        )
    

    def get_corpus(self, name):
        assert name == self.name
        return self
    

    def get_empirical_mutation_rate(self, use_weight=True, include_subclonal=True):

        # returns the ln mutation rate for each locus in the first sample
        kw=dict(use_weight = use_weight, include_subclonal=include_subclonal)
        mutation_rate = self.samples[0].get_empirical_mutation_rate(self.locus_dim, **kw)

        # loop through the rest of the samples and add the mutation rate using logsumexp
        for i in trange(1, len(self), desc = 'Piling up mutations', ncols=100):
            mutation_rate = mutation_rate + self.samples[i].get_empirical_mutation_rate(self.locus_dim, **kw)
        
        return mutation_rate
    

    def get_features_df(self):
        return DataFrame(
            {
                feature_name : feature['values']
                for feature_name, feature in self.features.items()
                if not feature['type'] == 'cardinality'
            }
        )
    
    def get_cardinality_features_df(self):
        return DataFrame(
            {
                feature_name : feature['values']
                for feature_name, feature in self.features.items()
                if feature['type'] == 'cardinality'
            }
        )



class MetaCorpus(CorpusMixin):
    
    def __init__(self, *corpuses):
        
        assert len(corpuses) > 1, 'If only one corpus, use that directly'
        assert len(set([corpus.name for corpus in corpuses])) == len(corpuses), \
            'All corpuses must have unique names.'
        assert all([np.all(corpuses[0].feature_names == corpus.feature_names) for corpus in corpuses])
        assert all([corpuses[0].shape == corpus.shape for corpus in corpuses])

        self._corpuses = corpuses
        self.idx_starts = np.cumsum(
            [0] + [len(corpus) for corpus in self.corpuses]
        )

    @property
    def observation_class(self):
        return self.corpuses[0].observation_class
    

    @property
    def corpuses(self):
        return self._corpuses

    @property
    def num_mutations(self):
        return sum([corpus.num_mutations for corpus in self.corpuses])

    @property
    def shape(self):
        return self.corpuses[0].shape
    
    @property
    def feature_dim(self):
        return self.corpuses[0].feature_dim
    
    @property
    def cardinality_features_dim(self):
        return self.corpuses[0].strand_dim
    
    @property
    def cardinalities_dim(self):
        return self.corpuses[0].cardinality_dim
    
    @property
    def locus_dim(self):
        return self.corpuses[0].locus_dim
    
    @property
    def corpus_names(self):
        return [corpus.name for corpus in self.corpuses]


    def __len__(self):
        return sum([len(corpus) for corpus in self.corpuses])


    def _get_corpus_idx(self, idx):

        for corpus_idx, (idx_start, next_idx) in enumerate(zip(self.idx_starts, self.idx_starts[1:])):
            
            if idx >= idx_start and idx < next_idx:
                return corpus_idx, idx - idx_start

        raise IndexError()


    def __getitem__(self, idx):
    
        corpus_idx, sample_idx = self._get_corpus_idx(idx)
        return self.corpuses[corpus_idx][sample_idx]


    def __iter__(self):

        for corpus in self.corpuses:
            for sample in corpus:
                yield sample


    def subset_samples(self, subset_idx):
        
        corpus_idxs = defaultdict(list)
        for s in subset_idx:
            corp, idx = self._get_corpus_idx(s)
            corpus_idxs[corp].append(idx)

        return MetaCorpus(
            *[
                self.corpuses[i].subset_samples(idxs)
                for i, idxs in corpus_idxs.items()
            ]
        )


    def subset_loci(self, loci):
        
        return MetaCorpus(*[
            corpus.subset_loci(loci) for corpus in self.corpuses
        ])
    

    def get_corpus(self, name):
        
        for corpus in self.corpuses:
            if corpus.name == name:
                return corpus
            
        raise KeyError(f'Corpus {name} does not exist.')

    @property
    def shared_correlates(self):
        return False
    
    @property
    def context_frequencies(self):
        return self.corpuses[0].context_frequencies

    @property
    def feature_names(self):
        return self.corpuses[0].feature_names



def train_test_split(corpus, seed = 0, train_size = 0.7,
                     by_locus = False):

    randomstate = np.random.RandomState(seed)

    N = len(corpus) if not by_locus else corpus.locus_dim
    subset_fn = corpus.subset_samples if not by_locus else corpus.subset_loci

    train_idx = sorted(randomstate.choice(
                N, size = int(train_size * N), replace=False
            ))
    
    test_idx = sorted(list( set(range(N)).difference(set(train_idx)) )) 

    return subset_fn(train_idx), subset_fn(test_idx)
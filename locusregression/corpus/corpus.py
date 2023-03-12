import numpy as np
from abc import ABC, abstractmethod
from collections import defaultdict
import h5py as h5
import logging
logger = logging.getLogger('Corpus')

class CorpusMixin(ABC):

    def __init__(self,*,
        name,
        samples,
        feature_names,
        X_matrix,
        trinuc_distributions,
        shared_exposures,
    ):
        self.name = name
        self.samples = samples
        self.feature_names = feature_names
        self.X_matrix = X_matrix
        self.trinuc_distributions = trinuc_distributions
        self._shared_exposures = shared_exposures

        if self._shared_exposures:
            self._exposures = samples[0]['exposures']

    @property
    def shape(self):
        return self.X_matrix.shape

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
    def get_corpus(self, name):
        raise NotImplementedError()



class SampleLoader:

    def __init__(self, filename, subset_idx = None):
        
        self.filename = filename
        
        if subset_idx is None:
            
            with h5.File(self.filename, 'r') as f:
                n_samples = len(f['samples'].keys())

            subset_idx = list(range(n_samples))
        
        self.subset_idx = subset_idx


    def __len__(self):
        return len(self.subset_idx)


    def _read_item(self, h5, i):
        #logger.debug('Streaming from disk cache.')
        return {k : h5[f'samples/{i}/{k}'][...] for k in h5[f'samples/{i}'].keys()}

    def __iter__(self):

        with h5.File(self.filename, 'r') as f:
            for i in self.subset_idx:
                yield self._read_item(f,i)
        

    def __getitem__(self, idx):
        
        idx = self.subset_idx[idx]

        with h5.File(self.filename, 'r') as f:
            return self._read_item(f, idx)


    def subset(self, idx_list):
        return SampleLoader(self.filename, [self.subset_idx[i] for i in idx_list])


class InMemorySamples(list):

    def subset(self, idx_list):
        return InMemorySamples([self[i] for i in idx_list])


def save_corpus(corpus, filename):

    with h5.File(filename, 'w') as f:
        
        data_group = f.create_group('data')
        data_group.attrs['shared_exposures'] = corpus.shared_exposures
        data_group.attrs['name'] = corpus.name

        data_group.create_dataset('trinuc_distributions', data = corpus.trinuc_distributions)
        data_group.create_dataset('X_matrix', data = corpus.X_matrix)
        data_group['X_matrix'].attrs['feature_names'] = corpus.feature_names

        samples_group = f.create_group('samples')

        for i, sample in enumerate(corpus.samples):
            for k, v in sample.items():
                samples_group.create_dataset(f'{i}/{k}', data = v)



def load_corpus(filename):

    with h5.File(filename, 'r') as f:

        is_shared = f['data'].attrs['shared_exposures']

        return Corpus(
            trinuc_distributions = f['data/trinuc_distributions'][...],
            X_matrix = f['data/X_matrix'][...],
            feature_names = f['data/X_matrix'].attrs['feature_names'],
            samples = InMemorySamples([
                {k : f[f'samples/{i}/{k}'][...] for k in f[f'samples/{i}'].keys()}
                for i in range(len(f['samples'].keys()))
            ]),
            shared_exposures=is_shared,
            name = f['data'].attrs['name']
        )


def stream_corpus(filename):

    with h5.File(filename, 'r') as f:

        is_shared = f['data'].attrs['shared_exposures']

        return Corpus(
            trinuc_distributions = f['data/trinuc_distributions'][...],
            X_matrix = f['data/X_matrix'][...],
            feature_names = f['data/X_matrix'].attrs['feature_names'],
            samples = SampleLoader(filename),
            shared_exposures=is_shared,
            name = f['data'].attrs['name']
        )


def train_test_split(corpus, seed = 0, train_size = 0.7):

    randomstate = np.random.RandomState(seed)

    train_idx = randomstate.choice(
        len(corpus), 
        size = int(train_size * len(corpus)),
        replace=False
    )

    test_idx = list( set( range(len(corpus)) ).difference(set(train_idx)) )

    return corpus.subset_samples(train_idx), corpus.subset_samples(test_idx)


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Corpus(CorpusMixin):

    def __len__(self):
        return len(self.samples)


    def __getitem__(self, idx):
        
        return dotdict({
            **self.samples[idx],
            'corpus_name' : self.name
        })


    def __iter__(self):

        for i in range(len(self)):
            yield self[i]

    @property
    def corpuses(self):
        return [self]

    def subset_samples(self, subset_idx):

        return Corpus(
            samples = self.samples.subset(subset_idx),
            feature_names = self.feature_names,
            X_matrix = self.X_matrix,
            trinuc_distributions = self.trinuc_distributions,
            shared_exposures = self.shared_exposures,
            name = self.name,
        )


    def subset_loci(self, loci):

        subsample_lookup = dict(zip(loci, np.arange(len(loci)).astype(int)))
        
        bool_array = np.zeros(self.X_matrix.shape[1]).astype(bool)
        bool_array[loci] = True

        total_mutations = 0
        new_samples = []
        for sample in self.samples:
                        
            mask = bool_array[sample['locus']]

            new_sample = {
                'mutation' : sample['mutation'][mask],
                'context' : sample['context'][mask],
                'count' : sample['count'][mask],
                'locus' : np.array([subsample_lookup[locus] for locus in sample['locus'][mask]]).astype(int), 
                'exposures' : sample['exposures'][:,loci]
            }

            total_mutations += new_sample['count'].sum()

            new_samples.append(new_sample)
        
        return Corpus(
            samples = InMemorySamples(new_samples),
            **{
                'X_matrix' : self.X_matrix[:,loci],
                'trinuc_distributions' : self.trinuc_distributions[:,loci],
                'feature_names' : self.feature_names,
            },
            shared_exposures = self.shared_exposures,
            name = self.name,
        )
    
    def get_corpus(self, name):
        assert name == self.name
        return self


class MetaCorpus(CorpusMixin):
    
    def __init__(self, *corpuses):
        
        assert len(corpuses) > 1, 'If only one corpus, use that directly'
        assert len(set([corpus.name for corpus in corpuses])) == len(corpuses), \
            'All corpuses must have unique names.'
        assert all([np.all(corpuses[0].feature_names == corpus.feature_names) for corpus in corpuses])
        assert all([corpuses[0].X_matrix.shape == corpus.X_matrix.shape for corpus in corpuses])

        self._corpuses = corpuses
        self.idx_starts = np.cumsum(
            [0] + [len(corpus) for corpus in self.corpuses]
        )

    @property
    def corpuses(self):
        return self._corpuses

    @property
    def shape(self):
        return self.corpuses[0].X_matrix.shape


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
    def trinuc_distributions(self):
        return self.corpuses[0].trinuc_distributions

    @property
    def feature_names(self):
        return self.corpuses[0].feature_names

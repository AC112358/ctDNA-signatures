from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
import h5py as h5


@dataclass
class Sample:
    """
    Represents a sample in the corpus.

    Attributes:
        attribute (np.ndarray): Array of attributes.
        mutation (np.ndarray): Array of mutations.
        context (np.ndarray): Array of contexts.
        locus (np.ndarray): Array of loci.
        name (str): Name of the sample.
        chrom (np.ndarray, optional): Array of chromosomes. Defaults to None.
        pos (np.ndarray, optional): Array of positions. Defaults to None.
        weight (np.ndarray, optional): Array of weights. Defaults to None.
    """

    attribute : np.ndarray
    mutation : np.ndarray
    context : np.ndarray
    locus : np.ndarray
    cardinality : np.ndarray
    name : str
    corpus_name : str
    corpus_name : str = None
    chrom : np.ndarray = None
    pos : np.ndarray = None
    weight : np.ndarray = None
    


    type_map = {
        'attribute' : np.uint8,
        'mutation' : np.uint8,
        'context' : np.uint8,
        'locus' : np.uint32,
        'cardinality' : np.uint8,
        'chrom' : 'S',
        'pos' : np.uint32,
        'weight' : np.float32,
        'name' : 'S',
        'corpus_name' : 'S',
    }

    data_attrs = ['attribute','mutation','context','locus','cardinality','weight','chrom','pos']
    required = ['attribute', 'mutation','context','locus','cardinality']
    mutation_attrs = ['attribute','mutation','context','locus','cardinality','weight','chrom','pos']

    N_CARDINALITY=1
    N_CONTEXTS=256
    N_MUTATIONS=3
    N_ATTRIBUTES=1

    def __len__(self):
        return len(self.locus)
    
    def __getitem__(self, i):
        return {
            k : getattr(self, k)[i]
            for k in self.mutation_attrs if hasattr(self, k)
        }
    
    @property
    def n_mutations(self):
        return sum(self.weight)

            
    def create_h5_dataset(self, h5_object, dataset_name):
        """
        Creates an HDF5 dataset from the sample data.

        Args:
            h5_object (h5py.File): The HDF5 file object.
            dataset_name (str): The name of the dataset.
        """
        for attr in self.data_attrs:
            h5_object.create_dataset(f'{dataset_name}/{attr}', data = self.__getattribute__(attr))

        h5_object[dataset_name].attrs['name'] = self.name

    
    @classmethod
    def read_h5_dataset(cls, h5_object, dataset_name, read_optional = True):
        """
        Reads an HDF5 dataset and returns a Sample object.

        Args:
            h5_object (h5py.File): The HDF5 file object.
            dataset_name (str): The name of the dataset.
            read_optional (bool, optional): Whether to read optional attributes. Defaults to True.

        Returns:
            Sample: The Sample object.
        """
        read_attrs = cls.data_attrs if read_optional else cls.required

        return cls(
            **{
                attr : h5_object[f'{dataset_name}/{attr}'][...]
                for attr in read_attrs
            },
            name = h5_object[dataset_name].attrs['name'],
        )
        
    
    def get_empirical_mutation_rate(self, n_loci, 
                                    use_weight = True, 
                                    include_subclonal=True):
        """
        Calculates the empirical mutation rate.

        Args:
            use_weight (bool, optional): Whether to use weights. Defaults to True.

        Returns:
            scipy.sparse.dok_matrix: The empirical mutation rate.
        """

        mutations = np.zeros((self.N_CARDINALITY, self.N_CONTEXTS, n_loci), dtype = float if use_weight else np.uint32)

        for mutation in self:
            if include_subclonal or mutation['weight'] >= 0.5:
                mutations[mutation['cardinality'], mutation['context'], mutation['locus']] += mutation['weight'] if use_weight else 1
        
        return mutations


    def plot_factorized(self, context_dist, mutation_dist, attribute_dist, ax=None, **kwargs):
        """
        Plots the data in from factorized form: (context, mutation, attribute).

        Args:
            ax: The matplotlib axes object.
            **kwargs: Additional keyword arguments for plotting.
        """
        pass

    
    def plot(self):
        pass

    
    def __getitem__(self, i):
        return {
            k : getattr(self, k)[i]
            for k in self.mutation_attrs if hasattr(self, k)
        }




class SampleLoader:
    """
    A class for lazy loading collections of samples from an h5 file.

    Args:
        filename (str): The path to the file containing the samples.
        subset_idx (list, optional): A list of indices specifying the subset of samples to load. 
            If not provided, all samples will be loaded.

    Attributes:
        filename (str): The path to the file containing the samples.
        subset_idx (list): A list of indices specifying the subset of samples to load.

    Methods:
        __len__(): Returns the number of samples in the loader.
        __iter__(): Returns an iterator over the samples.
        __getitem__(idx): Returns the sample at the specified index.
        subset(idx_list): Returns a new SampleLoader object with a subset of samples.
    """

    def __init__(self, filename, subset_idx=None, observation_class=Sample):
        self.filename = filename
        self.observation_class=observation_class

        if subset_idx is None:
            with h5.File(self.filename, 'r') as f:
                n_samples = len(f['samples'].keys())

            subset_idx = list(range(n_samples))

        with h5.File(self.filename, 'r') as f:
            all_sample_names = list(f['samples'].keys())
        
        self.sample_names = [all_sample_names[i] for i in subset_idx]
        self.subset_idx = subset_idx
        

    def __len__(self):
        return len(self.subset_idx)

    def _read_item(self, h5, idx):
        return self.observation_class.read_h5_dataset(h5, f'samples/{self.sample_names[idx]}')

    def __iter__(self):
        with h5.File(self.filename, 'r') as f:
            for i in range(len(self)):
                yield self._read_item(f, i)

    def __getitem__(self, idx):
        with h5.File(self.filename, 'r') as f:
            return self._read_item(f, idx)

    def subset(self, idx_list):
        return SampleLoader(self.filename, [self.subset_idx[i] for i in idx_list])


class InMemorySamples(list):

    def subset(self, idx_list):
        return InMemorySamples([self[i] for i in idx_list])
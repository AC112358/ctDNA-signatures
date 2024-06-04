import numpy as np
from ..corpus.sbs.observation_config import SBSSample, CONTEXT_IDX, MUTATIONS_IDX
from ..corpus.length.observation_config import LENGTH_BINS
from ..corpus import BED12Record, Corpus, InMemorySamples
import json
import os
import tqdm
from joblib import Parallel, delayed
from functools import partial
import pandas as pd
from itertools import product

TRANSITION_MATRIX = np.array([
    [0.99, 0.005, 0.005],
    [0.0, 0.97, 0.03],
    [0.015, 0.015, 0.97]
])

BETA_MATRIX = np.array([
    [1,-0.5,-0.5],
    [-0.5,1,-0.5],
    [-0.5,-0.5,1],
    [1,-1,1],
])

TRINUC_PRIORS = np.array([
    np.ones(32) * 5.,
    np.ones(32) * 20,
    np.ones(32) * 10,
])

SIGNAL_MEANS = np.array([1.,1.,1.])
SIGNAL_STD = np.array([0.3, 0.25, 0.5])

complement = {'A' : 'T','T' : 'A','G' : 'C','C' : 'G'}

def revcomp(seq):
    return ''.join(reversed([complement[nuc] for nuc in seq]))


with open(os.path.join(os.path.dirname(__file__), 'cosmic.json'),'r') as f:
    COSMIC_SIGS = json.load(f)

# IN5P_SIGS = pd.read_csv(os.path.join(os.path.dirname(__file__),'MSK_allreads_W_in5p_NMF_10_20_17.csv'), index_col=0) # without normalization
IN5P_SIGS = pd.read_csv(os.path.join(os.path.dirname(__file__),'W_in5p_catalog_in_v2.csv'), index_col=0) # re-unnormalized
IN5P_SIGS = IN5P_SIGS.to_dict()
OUT5P_SIGS = pd.read_csv(os.path.join(os.path.dirname(__file__),'W_out_sparse_selected.csv'), index_col=0) # re-unnormalized
OUT5P_SIGS = OUT5P_SIGS.to_dict()
LEN_SIGS = pd.read_csv(os.path.join(os.path.dirname(__file__),'W_length_k10.csv'), index_col=0) # re-unnormalized
LEN_SIGS = LEN_SIGS.to_dict()


class SimulatedCorpus:

    @staticmethod
    def cosmic_sig_to_matrix(cosmic_sig):
        
        sigmatrix = np.zeros((len(CONTEXT_IDX),3))

        for key, p in cosmic_sig.items():
            context = key[0] + key[2] + key[6]
            mutation = key[4]

            sigmatrix[CONTEXT_IDX[context], MUTATIONS_IDX[context][mutation]] += p

        return sigmatrix

    @staticmethod
    def motif_sig_to_matrix(motif_sig):
        CONTEXTS = sorted(
            map(lambda x : ''.join(x), product('ATCG','ATCG','ATCG', 'ATCG')), 
            key = lambda x : (x[0], x[1], x[2], x[3])
            )
        
        CONTEXT_IDX = dict(zip(CONTEXTS, range(len(CONTEXTS))))
        sigmatrix = np.zeros((len(CONTEXT_IDX),1))
    
        for context, p in motif_sig.items():
            sigmatrix[CONTEXT_IDX[context], 0] = p
    
        return sigmatrix
    
    def len_sig_to_matrix(length_sig): # I haven't tested this yet, I have to understand why the matrix has context and mutation dimensions rather than context x mutation vs sample

        sigmatrix = np.zeros((1, len(LENGTH_BINS)))
        
        length_bin_names = [f"{min_length}-{max_length}" for min_length, max_length in LENGTH_BINS]
        mutation_idx = {name: idx for idx, name in enumerate(length_bin_names)}

        for length_bin, p in length_sig.items():
            sigmatrix[0, mutation_idx[length_bin]] = p

        return sigmatrix

    
    @staticmethod
    def create(
        corpus_name,
        seed = 0,
        n_cells = 100,
        log_mean_mutations = 5,
        log_std_mutations = 1.,
        pi_prior = 1.,
        n_loci = 1000,
        mutation_rate_noise = 0.1,
        state_transition_matrix = TRANSITION_MATRIX.copy(),
        beta_matrix = BETA_MATRIX.copy(),
        rate_function = None,
        trinucleotide_priors = TRINUC_PRIORS.copy(),
        signal_means = SIGNAL_MEANS.copy(),
        signal_std = SIGNAL_STD.copy(),
        exposures = None,*,
        cosmic_sigs,
        cardinality_effects,
    ):

        signatures = np.vstack([
            np.expand_dims(SimulatedCorpus.cosmic_sig_to_matrix(COSMIC_SIGS[sig_id]), 0)
            for sig_id in cosmic_sigs
        ])

        num_states = state_transition_matrix.shape[0]
        assert state_transition_matrix.shape == (num_states, num_states)

        num_signatures = signatures.shape[0]
        assert beta_matrix is None or beta_matrix.shape == (num_signatures, num_states)

        if beta_matrix is None:
            assert not rate_function is None

        #assert signatures.shape == (num_signatures, 32, 3)

        assert isinstance(signal_means, np.ndarray) and isinstance(signal_std, np.ndarray) and \
            signal_means.shape == (num_states,) and signal_std.shape == (num_states,)

        assert trinucleotide_priors.shape == (num_states, len(CONTEXT_IDX))

        #assert isinstance(pi_prior, (int, float))
        if exposures is None:
            exposures = np.ones((n_cells, n_loci))
        else:
            exposures = np.ones((1, n_loci))

        randomstate = np.random.RandomState(seed)

        states = SimulatedCorpus.get_genomic_states(randomstate, 
            n_loci=n_loci, transition_matrix=state_transition_matrix
        )

        cardinality_states = SimulatedCorpus.get_genomic_states(randomstate, 
            n_loci=n_loci, transition_matrix=state_transition_matrix
        )

        signals = SimulatedCorpus.get_signals(randomstate, state = states, 
            signal_means=signal_means, signal_std=signal_std
        )

        # 2 x C x L
        context_frequencies = np.vstack([
            randomstate.dirichlet(np.tile(trinucleotide_priors[state], 2))[None,:]
            for state in states
        ]).T.reshape(2, -1, n_loci)

        omega = signatures
        delta = signatures.sum(axis=-1)/context_frequencies.sum(axis=(0,2))

        cell_pi = randomstate.dirichlet(np.ones(num_signatures) * pi_prior, size = n_cells)
        cell_n_mutations = randomstate.lognormal(log_mean_mutations, log_std_mutations, size = n_cells).astype(int)

        locus_effects = SimulatedCorpus.get_locus_effects(
                    signals, 
                    exposures, 
                    beta_matrix = beta_matrix, 
                    rate_function = rate_function,
                    mutation_rate_noise = mutation_rate_noise,
                    random_state= randomstate,
                )
        
        psi_matrix = SimulatedCorpus._get_psi_matrix(
            exposures = exposures,
            context_frequencies = context_frequencies,
            delta = delta,
            locus_effects = locus_effects,
            cardinality_states = cardinality_states,
            cardinality_effects = cardinality_effects
        )

        samples = []
        for i, (pi, n_mutations) in tqdm.tqdm(enumerate(zip(cell_pi, cell_n_mutations)),
            total = len(cell_pi), ncols = 100, desc = 'Generating samples'):

            samples.append(
                SimulatedCorpus.simulate_sample(
                    randomstate = randomstate,
                    omega = omega,
                    pi = pi,
                    n_mutations = n_mutations,
                    exposures = exposures,
                    name = i,
                    psi_matrix = psi_matrix
                )
            )

        feature_dict = {
                    'feature' + str(k) : {
                        'type' : 'continuous',
                        'values' : signal,
                        'group' : 'all'
                    }
                    for k, signal in enumerate(signals)
                }
        
        card_map = ['-','.','+']
        feature_dict['cardinality'] = {
            'type' : 'cardinality',
            'values' : np.array([card_map[c] for c in cardinality_states]),
            'group' : 'all'
        }

        corpus = Corpus(
                type='SBS',
                name = corpus_name,
                samples = InMemorySamples(samples),
                context_frequencies = context_frequencies,
                features = feature_dict,
                regions = [
                    BED12Record(1, 0, 1000, str(i), 0, '+', 0, 1000, '0,0,0', 1, [1000], [0])
                    for i in range(n_loci)
                ]
        )

        generative_parameters = {
            'states' : states,
            'compositions' : cell_pi,
            'beta' : beta_matrix,
            'signatures' : signatures,
            'n_cells' : n_cells,
            'log_mean_mutations' : log_mean_mutations,
            'log_std_mutations' : log_std_mutations,
            'n_loci' : n_loci,
            'pi_prior' : pi_prior,
            'mutation_rate_noise' : mutation_rate_noise,
            'seed' : seed,
            'cosmic_sigs' : cosmic_sigs
        }

        return corpus, generative_parameters


    @staticmethod
    def get_genomic_states(randomstate,*, n_loci, transition_matrix):
        
        state = 0
        history = []
        for _ in range(n_loci):

            state = randomstate.choice(
                len(transition_matrix[state]),
                p = transition_matrix[state]
            )

            history.append(state)
            
        return np.array(history)


    @staticmethod
    def get_signals(randomstate,*, state, signal_means, signal_std):
        
        n_loci, n_features = len(state), len(signal_means)
        state_as_onehot = np.zeros(( n_loci, n_features ))
        state_as_onehot[ np.arange(n_loci), state] = np.ones(n_loci).astype(int)

        signals = signal_means[:,None]*state_as_onehot.T \
                + signal_std[:,None]*randomstate.randn(len(signal_means), len(state)) # S, L

        standard_signals = (signals - signals.mean(-1, keepdims = True))/signals.std(-1, keepdims = True)

        return standard_signals


    @staticmethod
    def get_locus_effects(signals, exposure,
                        beta_matrix = None, 
                        rate_function = None,
                        mutation_rate_noise = 0.5,*,
                        random_state):
        
        if rate_function is None:
            assert not beta_matrix is None
            rate_function = lambda x : beta_matrix @ x

        return np.exp(rate_function(signals) + mutation_rate_noise*random_state.randn(signals.shape[-1]))


    @staticmethod
    def _get_psi_matrix(*, exposures, context_frequencies, delta, locus_effects, cardinality_effects, cardinality_states):

        locus_effects = (np.log(locus_effects) + np.log(exposures)) # KxL + 1xL -> KxL

         # Kx1xCx1 + 1x2xCxL -> Kx2xCxL
        signature_effects = np.log(delta)[:,None,:,None] + np.log(context_frequencies)[None,:,:,:] + np.log(10000) \
                                + cardinality_states[None,None,None,:]*np.log(cardinality_effects)[:,None,None,None]*np.array([-1,1])[None,:,None,None]

        logits = locus_effects[:, None, None, :] + signature_effects # Kx2xCxL

        logdenom = np.logaddexp.reduce(logits, axis = (1,2,3), keepdims = True) # Kx1x1x1

        return np.exp(
            np.nan_to_num(
                logits - logdenom,
                nan = -np.inf
            )
        )
    

    @staticmethod
    def simulate_sample(
            randomstate = None,
            seed = None,*,
            omega,
            pi, 
            n_mutations, 
            exposures, 
            name,
            psi_matrix,
        ):

        if randomstate is None:
            randomstate = np.random.RandomState(seed)
        
        loci=[]; contexts=[]; mutations=[]; cardinalities=[]; z=[]
        # KxDxCxL --> KxL
        p_l = psi_matrix.sum(axis = (1,2))

        # KxDxCxL --> KxDxL
        p_card_given_l = psi_matrix.sum(axis = 2)/p_l[:,None,:]
        
        for i in range(n_mutations):

            sig = randomstate.choice(len(pi), p = pi)

            locus = randomstate.choice(p_l.shape[1], p = p_l[sig,:])
            cardinality = randomstate.choice([0,1], p = p_card_given_l[sig,:,locus])

            context_dist = psi_matrix[sig,cardinality,:,locus]/psi_matrix[sig,cardinality,:,locus].sum()
            
            context = randomstate.choice(len(CONTEXT_IDX), p = context_dist)

            mutation_dist = omega[sig, context]/omega[sig, context].sum()
            mutation = randomstate.choice(3, p = mutation_dist)

            loci.append(locus)
            contexts.append(context)
            mutations.append(mutation)
            cardinalities.append(cardinality)
            z.append(sig)


        elements = dict(
            mutation=np.array(mutations), 
            context=np.array(contexts), 
            locus=np.array(loci), 
            weight=np.ones_like(loci).astype(float),
            cardinality=np.array(cardinalities),
            pos = np.array(z),
            chrom = np.array(['chr1']*len(loci)),
            attribute = np.zeros_like(loci),
        )

        for k, v in elements.items():
            elements[k] = v.astype(SBSSample.type_map[k])

        return SBSSample(
                        **elements,
                        name = str(name),
                    )


    @staticmethod
    def from_model(
        model, 
        corpus_state,
        corpus,
        use_signatures = None,
        seed = 0,
        n_jobs = 1,
        ):

        randomstate = np.random.RandomState(seed)
        n_loci = corpus.locus_dim

        exposures = np.ones((1, n_loci))

        if not use_signatures is None:
            use_signatures = np.array([model.component_names.index(sig) for sig in use_signatures])
        else:
            use_signatures = np.arange(model.n_components)
            
        psi_matrix = np.exp(corpus_state.get_log_component_effect_rate(model.model_state, exposures))
        psi_matrix = psi_matrix[use_signatures]

        cell_gamma = np.array([
            model._predict_sample(sample, corpus_state)
            for sample in tqdm.tqdm(corpus, ncols = 100, desc = 'Calculating exposures')
        ])

        cell_gamma = cell_gamma[:,use_signatures]
        cell_pi = cell_gamma / cell_gamma.sum(axis = 1)[:,None]

        cell_n_mutations = cell_gamma.sum(axis = 1).astype(int)

        generate_sample_fn = partial(
            SimulatedCorpus.simulate_sample,
            omega=model.model_state._rho[use_signatures],
            exposures=exposures,
            psi_matrix=psi_matrix,
        )

        samples = Parallel(n_jobs=n_jobs, verbose = 10)(
            delayed(generate_sample_fn)(
                pi = cell_pi[sample_num],
                n_mutations = cell_n_mutations[sample_num],
                name = str(sample_num),
                seed = seed + sample_num,
            )
            for sample_num in range(len(cell_pi))
        )
        
        return Corpus(
                type='SBS',
                name = corpus.name,
                samples = InMemorySamples(samples),
                context_frequencies = corpus.context_frequencies,
                features = corpus.features,
                regions = corpus.regions
            )

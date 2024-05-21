import numpy as np
from scipy.spatial.distance import cdist
from pandas import DataFrame

def posterior_divergence(*,
        model_state, 
        sample, 
        corpus_state, 
        component_names, 
        n_iters=1000, 
        warmup=300,
):
    
    from ..model import get_sample_posterior
    
    def get_logit_loss(posterior_df):

        # the first two columns are CHR and POS
        log_q = posterior_df[posterior_df.columns[2:]].values[
                    np.arange(len(posterior_df)),
                    posterior_df.POS.values,
                ]

        return log_q
    
    z_posterior = get_sample_posterior(
        model_state=model_state,
        component_names=component_names,
        sample=sample,
        quiet=True,
        corpus_state=corpus_state,
        n_iters = n_iters,
        warmup = warmup,
    )

    return list(zip(
        [component_names[j] for j in z_posterior.POS.values], 
        get_logit_loss(z_posterior)
    ))


def signature_cosine_distance(model, simulation_parameters):
    
    sigs = np.vstack([
        model.signature(i, raw = True, normalization='global') for i in range(model.n_components)
    ])

    truth = simulation_parameters['signatures'].reshape(-1,96)

    cosine_matches = 1 - cdist(sigs, truth, metric='cosine')
    
    return cosine_matches.max(0).mean()



def coef_l1_distance(model, simulation_parameters):
    
    return cdist(model.model_state.beta_mu[:,:-1], 
            simulation_parameters['beta'], 
            metric='cityblock'
            ).min(0).mean()

import logging
from numpy.random import RandomState
from numpy import vstack, squeeze
from joblib import Parallel, delayed
import shap
import pandas as pd
logger = logging.getLogger(__name__)

def explain(
    signature,*,
    model,
    corpus,
    n_jobs = 1,
    chunk_size = None,
):
    
    def _calculate_shap_values(tree_explainer, chunk):
        return squeeze(
            tree_explainer.shap_values(
                chunk,
                check_additivity=False,
                approximate=False,
            )
        )
    
    component = model._get_signature_idx(signature)

    if not model.is_trained:
        logger.warn('This model was not trained to completion, results may be innaccurate')
    
    tree_model = model.model_state.rate_models[component]

    corpus_states = {
        name : state.clone_corpusstate(corpus.get_corpus(name))
        for name, state in model.corpus_states.items()
        if name in corpus.corpus_names
    }
    
    X_tild = model.model_state.feature_transformer\
        .transform(corpus_states)

    background_idx = RandomState(0)\
                        .choice(
                            len(X_tild), 
                            size = 1000, 
                            replace = False
                        )

    explainer = shap.TreeExplainer(
        tree_model,
        X_tild[background_idx],
    )

    shap_values = _calculate_shap_values(explainer, X_tild)

    feature_names = model.model_state.feature_transformer.feature_names_out

    return (
        shap_values,
        X_tild, 
        feature_names,
        model.model_state.feature_transformer.assemble_matrix(corpus_states)[feature_names].values,
    )

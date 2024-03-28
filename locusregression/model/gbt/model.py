from ..model import LocusRegressor
from ._gbt_modelstate import GBTModelState, GBTCorpusState
import locusregression.model.gbt._gbt_sstats as _gbt_sstats
import matplotlib.pyplot as plt


class GBTRegressor(LocusRegressor):

    MODEL_STATE = GBTModelState
    CORPUS_STATE = GBTCorpusState
    SSTATS = _gbt_sstats

    def __init__(self,
                 tree_learning_rate=0.1, 
                 max_depth = 5,
                 max_trees_per_iter = 25,
                 max_leaf_nodes = 31,
                 min_samples_leaf = 30,
                 max_features = 1.,
                 l2_regularization=1,
                 n_iter_no_change=2,
                 use_groups=True,
                 signature_reg = 0.,
                 cardinality_reg = 0.,
                  **kw, 
                ):
        super().__init__(**kw)
        self.tree_learning_rate = tree_learning_rate
        self.max_depth = max_depth
        self.use_groups = use_groups
        self.max_trees_per_iter = max_trees_per_iter
        self.l2_regularization=l2_regularization
        self.n_iter_no_change=n_iter_no_change
        self.max_leaf_nodes = max_leaf_nodes
        self.min_samples_leaf = min_samples_leaf
        self.max_features = max_features
        self.cardinality_reg = cardinality_reg
        self.signature_reg = signature_reg

    
    def _get_rate_model_parameters(self):
        return {
            'use_groups' : self.use_groups,
            'tree_learning_rate' : self.tree_learning_rate,
            'max_depth' : self.max_depth, 
            'max_trees_per_iter' : self.max_trees_per_iter,
            'n_iter_no_change' : self.n_iter_no_change,
            'l2_regularization': self.l2_regularization,
            'max_leaf_nodes' : self.max_leaf_nodes,
            'min_samples_leaf' : self.min_samples_leaf,
            'max_features' : self.max_features,
            'cardinality_reg' : self.cardinality_reg,
            'signature_reg' : self.signature_reg,
        }
    
    @classmethod
    def sample_params(cls, trial):
        return {
            'l2_regularization': trial.suggest_categorical('l2_regularization', [0.,1e-2,1e-1, 1, 10, 100, 1000]),
        }
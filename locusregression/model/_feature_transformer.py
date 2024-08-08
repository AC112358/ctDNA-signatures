import pandas as pd
from collections import defaultdict
from sklearn.preprocessing import OneHotEncoder, PowerTransformer, MinMaxScaler, \
    QuantileTransformer, StandardScaler, RobustScaler, LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.base import clone, BaseEstimator
from sklearn import set_config
import logging
from numpy import array, vstack
logger = logging.getLogger(' LocusRegressor')
set_config(enable_metadata_routing=True)


def _assemble_matrix(feature_names, features):
    '''
    For a single corpus, assemble a matrix of features.
    '''
    return pd.DataFrame(
            {feature_name : features[feature_name]['values'] 
             for feature_name in feature_names
            }
        )


class FeatureTransformer:

    def __init__(self, categorical_encoder = OneHotEncoder(sparse_output=False, drop='first')):
        self.categorical_encoder = clone(categorical_encoder)


    def list_feature_groups(self):
        
        feature_groups_dict_ = defaultdict(list)
        for feature_name, groups in zip(
            self.feature_names_, self.groups_
        ):
            for group in groups:
                try:
                    feature_groups_dict_[group].append(
                        self.feature_names_out.index(feature_name)
                    )
                except ValueError:
                    feature_groups_dict_[group].extend(
                        [idx for idx, name in enumerate(self.feature_names_out) if name.startswith(feature_name + '_')]
                    )

        return list(set(map(tuple,feature_groups_dict_.values())))


    @property
    def feature_names_out(self):
        return list(self.transformer_.get_feature_names_out())
    

    def list_categorical_features(self):
        #encoder_name = type(self.categorical_encoder).__name__.lower()
        _slice = self.transformer_.output_indices_['categorical']
        
        return list(range(_slice.start, _slice.stop))
    

    @property
    def transformer_(self):
        return next(iter(self.transformers_.values()))
    

    def get_base_transformer(self, corpus_states):

        example_features = next(iter(corpus_states.values())).features
        
        self.feature_names_ = list([f for f in example_features.keys() if example_features[f]['type'] != 'cardinality'])
        self.feature_types_ = [example_features[feature]['type'] for feature in self.feature_names_]
        self.groups_ = [example_features[feature]['group'].split(',') for feature in self.feature_names_]
        
        self.feature_type_dict_ = defaultdict(list)        
        for idx, feature_type in enumerate(self.feature_types_):
            self.feature_type_dict_[feature_type].append(idx)

        # Set the categories for the categorical encoder so that it uses the 
        # same categories across all datasets
        example_matrix = _assemble_matrix(self.feature_names_, example_features)
        categorical_features = example_matrix.iloc[:, self.feature_type_dict_['categorical']]

        cat_encoder = clone(self.categorical_encoder).set_params(
            categories=[
                list(categorical_features[feature].unique()) 
                for feature in categorical_features.columns
            ])

        return ColumnTransformer([
            ('power', PowerTransformer(), self.feature_type_dict_['power']),
            ('minmax', MinMaxScaler(), self.feature_type_dict_['minmax']),
            ('quantile', QuantileTransformer(output_distribution='uniform'), self.feature_type_dict_['quantile']),
            ('standardize', StandardScaler(), self.feature_type_dict_['standardize']),
            ('robust', RobustScaler(), self.feature_type_dict_['robust']),
            ('categorical', cat_encoder, self.feature_type_dict_['categorical']),],
            remainder='passthrough',
            verbose_feature_names_out=False,
        )


    def fit(self, corpus_states):
        
        self.transformers_ = {}
        base_transformer = self.get_base_transformer(corpus_states)

        for corpus_name, corpus_state in corpus_states.items():
            matrix = _assemble_matrix(self.feature_names_, corpus_state.features)
            self.transformers_[corpus_name] = clone(base_transformer).fit(matrix)

        # outro messages
        if len(self.feature_type_dict_['categorical']) > 0:
            logger.info(f'Found categorical features: {", ".join(self.feature_names_out[i] for i in self.list_categorical_features())}')
            
        for feature_group in self.list_feature_groups():
            if len(feature_group) > 0:
                logger.info(f'Found feature group: {", ".join(self.feature_names_out[i] for i in feature_group)}')

        return self


    def transform(self, corpus_states):

        for corpus_state in corpus_states.values():
            assert all([f in corpus_state.feature_names for f in self.feature_names_])

        if len(corpus_states) == 1:
            corpus_name, state = next(iter(corpus_states.items()))
            return self.transformers_[corpus_name].transform(
                _assemble_matrix(self.feature_names_, state.features)
            )
        else:
            #raise NotImplementedError("This path should not be used...")
            return vstack([
                self.transformers_[corpus_name].transform(
                    _assemble_matrix(self.feature_names_, corpus_state.features)
                )
                for corpus_name, corpus_state in corpus_states.items()
            ])
    


class CardinalityTransformer(BaseEstimator):

    def fit(self, corpus_states):
        example_features = next(iter(corpus_states.values())).features
        self.feature_names_ = list([f for f in example_features.keys() if example_features[f]['type'] == 'cardinality'])

        logger.info(
            f'Found cardinality features: {", ".join(self.feature_names_)}'
        )
        
        # theres actually no fitting to be done.
        #self.transformer_ = LabelEncoder(classes=['-','.','+']).fit(['-','.','+'])
        
        return self
    

    def transform(self, corpus_states):

        for corpus_state in corpus_states.values():
            assert all([f in corpus_state.feature_names for f in self.feature_names_])

        assert len(corpus_states) == 1

    
        _, state = next(iter(corpus_states.items()))
        matrix = _assemble_matrix(self.feature_names_, state.features)
        #matrix = matrix.values
        
        #transformed_matrix = array([
        #    self.transformer_.transform(matrix[:,col]) - 1
        #    for col in range(matrix.shape[1])
        #]).T
        conversion_dict={'+' : 1, '-' : -1, '.' : 0}
        transformed_matrix=matrix.map(lambda x : conversion_dict[x]).values

        return transformed_matrix
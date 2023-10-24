from ..model import LocusRegressor
from ._gbt_modelstate import GBTModelState, GBTCorpusState
import locusregression.model.gbt._gbt_sstats as _gbt_sstats

class GBTRegressor(LocusRegressor):

    MODEL_STATE = GBTModelState
    CORPUS_STATE = GBTCorpusState
    SSTATS = _gbt_sstats

    @classmethod
    def sample_params(cls, randomstate):
        return dict(
            seed = randomstate.randint(0, 100000000),
            tau = randomstate.choice([1,4,8, 16, 48, 128]),
            kappa = randomstate.choice([0.5, 0.6, 0.7, 0.8]),
        )

import numpy as np
from scipy.special import psi, gammaln, polygamma, gamma
import logging
logger = logging.getLogger(__name__)


def dirichlet_multinomial_logprob(z, alpha):

    z, num = np.unique(z, return_counts=True)

    n_z = np.zeros_like(alpha)
    n_z[z] = num

    n = sum(n_z)

    alpha_bar = sum(alpha)

    return gammaln(alpha_bar) + gammaln(n+1) - gammaln(n+alpha_bar) + \
            np.sum(
                gammaln(n_z + alpha) - gammaln(alpha) - gammaln(n_z + 1)
            )


def log_dirichlet_expectation(alpha):
    
    if len(alpha.shape) == 1:
        return psi(alpha) - psi(np.sum(alpha))
    else:
        return psi(alpha) - psi(np.sum(alpha, axis = -1, keepdims=True))


def dirichlet_bound(alpha, gamma):

    logE_gamma = log_dirichlet_expectation(gamma)

    alpha = np.expand_dims(alpha, 0)
    
    return gammaln(np.sum(alpha)) - gammaln(np.sum(gamma, axis = -1)) + \
        np.sum(
             gammaln(gamma) - gammaln(alpha) + (alpha - gamma)*logE_gamma,
             axis = -1
        )


def _dir_prior_update_step(prior, N, logphat):
    
    gradf = N * (psi(np.sum(prior)) - psi(prior) + logphat)

    c = N * polygamma(1, np.sum(prior))
    q = -N * polygamma(1, prior)

    b = np.sum(gradf / q) / (1 / c + np.sum(1 / q))

    dprior = -(gradf - b) / q

    return dprior + prior



def update_dir_prior(prior, N, logphat):

    def _check_prior(prior):
        return np.all(prior > 0) and np.all(np.isfinite(prior))

    old_prior = prior.copy()
    prior = _dir_prior_update_step(prior, N, logphat)

    if not _check_prior(prior):
        logger.warning(' Prior update failed. Retrying with different initialization.')
        prior = np.exp(logphat)

        prior = _dir_prior_update_step(prior, N, logphat)
        if not _check_prior(prior):
            logger.warning(' Prior update failed, reverting to old prior.')
            prior = old_prior
    
    return prior


def update_tau(mu, nu):
    return np.sqrt(2*np.pi) * np.sum(mu**2 + nu**2, axis = -1)


def update_alpha(alpha, gamma):
    
    N = gamma.shape[0]
    log_phat = log_dirichlet_expectation(gamma).mean(-2)
    return update_dir_prior(alpha, N, log_phat)
"""
Module contains all functions required to perform a single property
maximum loglike optimization.

Inspired/imported from:
* https://github.com/THGLab/X-EISD/blob/master/eisd/scorers.py
* https://github.com/Oufan75/X-EISD/blob/master/eisd/scorers.py
"""
import numpy as np


def calc_opt_params(beta, exp, exp_sig, sig):
    opt_params = np.zeros(beta.shape)
    if not np.any(exp_sig==0):
        ratio = (sig ** 2.0) / (exp_sig ** 2.0)
        opt_params = (ratio * (exp - beta)) / (1.0 + ratio)
    return opt_params


def normal_loglike(x, mu, sig, gamma=1.0):
    # Allow for cases where one/more sig are zero
    # all assigned to zero for now
    # modified by Oufan Zhang @Oufan75
    logp = np.zeros(x.shape)
    if not np.any(sig==0):
        exp_val = -gamma * ((x - mu) ** 2.0) / (2.0 * (sig ** 2.0))
        pre_exp = 1.0 / (np.sqrt(2.0 * np.pi * (sig ** 2.0)))
        logp = np.log(pre_exp * np.exp(exp_val))
    return logp


def calc_score(beta, exp, exp_sig, sig, opt_params, gamma=1.0):
    f_q = normal_loglike(opt_params, 0, sig, gamma)
    err = exp - opt_params - beta
    f_err = normal_loglike(err, 0, exp_sig, gamma)
    f = f_q + f_err
    f_comps = [f_q, f_err]
    return f, f_comps


def get_gamma(exp_saxs, num_res, qmax, qmin):
    # Generalizable way to calculate gamma parameter
    # 
    N_q = len(exp_saxs)
    D_max = num_res
    N_s = D_max * (qmax - qmin) / np.pi
    
    return N_s / N_q


def saxs_optimization_ensemble(
    q_max,
    q_min,
    exp_data,
    bc_data,
    indicies,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    #TODO: complete SAXS module
    return total_score, error, rmse, bc_saxs
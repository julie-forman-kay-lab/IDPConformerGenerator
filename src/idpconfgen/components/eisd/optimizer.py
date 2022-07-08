"""
Inspired/imported from https://github.com/THGLab/X-EISD/blob/master/eisd/optimizer.py
"""
import numpy as np

from idpconfgen.components.eisd import modes, eisd_run_all
from idpconfgen.components.eisd.scorers import *


def core_eisd(
    exp_data,
    bc_data,
    ens_size,
    epochs,
    pre_idx=None,
    mode=eisd_run_all,
    beta=0.1,
    opt_type=None,
    ):
    """
    Main function to calculate scores and potentially optimize ensembles
    given experimental and back-calculated data.

    Parameters
    ----------
    exp_data : dict
        Dictionary to the paths of experimental data files for each
        module being the key.
    
    bc_data : dict
        Dictionary to the paths of back-calculated data files for each
        module being the key.
        
    ens_size : int
        Number of conformers in the ensemble.
    
    epochs : int
        Number of times to run main optimization.
    
    pre_idx : ndarray, optional
        Shape (number of ensembles, size of ensemble). Fastest way
        
    mode : str
        Which EISD scoring modes to run.
        Defaults to `all`.
    
    beta : float
        Hyperparameter for monte-carlo type optimizations.
        Related to temperature. Defaults to 0.1.

    opt_type : str
        Optimization type can be `mc` for Metropolis Monte Carlo
        or `max` for score optimization method. If None, unoptimized
        results will be returned.
    """
    # get size of pool
    pool_size = bc_data[next(iter(bc_data))].data.shape[0]
    
    # switch the property
    flags = modes(mode)
    
    final_results = []
    final_indices = []
    final_best_jcoups = []
    EMPTY = [0, 0, 0]
    
    if pre_idx:
        epochs = pre_idx.shape[0]
        flags = modes(eisd_run_all)
        ens_size = pre_idx.shape[1]
    
    for i in range(epochs):
        # random selection without replacement
        if pre_idx is None:
            indices = np.random.choice(np.arange(pool_size), ens_size, replace=False)
        else:
            indices = pre_idx[i]

        if flags[saxs_name]:
            rmse_saxs, total_score_saxs, bc_saxs, _ = saxs_optimization_ensemble(exp_data, bc_data, indices)
        else:
            rmse_saxs, total_score_saxs, bc_saxs = EMPTY

        if flags[cs_name]:
            rmse_cs, total_score_cs, bc_cs, _ = cs_optimization_ensemble(exp_data, bc_data, indices)
        else:
            cs_num = 1.
            rmse_cs, total_score_cs, bc_cs = EMPTY

        if flags[fret_name]:
            rmse_fret, total_score_fret, bc_fret, _ = fret_optimization_ensemble(exp_data, bc_data, indices)
        else:
            rmse_fret, total_score_fret, bc_fret = EMPTY
        #TODO: double check with Oufan to see if we should use jc_ensemble()
        if flags[jc_name]:
            rmse_jc, total_score_jc, best_jcoups, _, bc_alpha_vals_jc = jc_optmization_ensemble(exp_data, bc_data, indices)
        else:
            jc_num = 1.
            rmse_jc, total_score_jc, bc_alpha_vals_jc, best_jcoups = [0, 0, 0, [0]]

        if flags[noe_name]:
            noe_num = exp_data[noe_name].data.shape[0]
            rmse_noe, total_score_noe, bc_dist_vals_noe, _ = noe_optimization_ensemble(exp_data, bc_data, indices)
        else:
            noe_num = 1.
            rmse_noe, total_score_noe, bc_dist_vals_noe = EMPTY

        if flags[pre_name]:
            rmse_pre, total_score_pre, bc_dist_vals_pre, _  = pre_optimization_ensemble(exp_data, bc_data, indices)
        else:
            pre_num = 1.
            rmse_pre, total_score_pre, bc_dist_vals_pre = EMPTY

        if flags[rdc_name]:
            rmse_rdc, total_score_rdc, bc_rdc, _  = rdc_optimization_ensemble(exp_data, bc_data, indices)
        else:
            rmse_rdc, total_score_rdc, bc_rdc = EMPTY

        if flags[rh_name]:
            rmse_rh, total_score_rh, bc_rh, _  = rh_optimization_ensemble(exp_data, bc_data, indices)
        else:
            rmse_rh, total_score_rh, bc_rh = EMPTY

        if opt_type:
            max_score_SAXS, opt_rmse_SAXS, max_score_CS, opt_rmse_CS, \
            max_score_FRET, opt_rmse_FRET, max_score_JC, opt_rmse_JC, \
            max_score_NOEs, opt_rmse_NOE, max_score_PREs, opt_rmse_PRE, \
            max_score_RDCs, opt_rmse_RDC, max_score_RH, opt_rmse_RH, \
            accepted, new_indices, best_jcoups = maximize_score(
                exp_data, bc_data, ens_size, indices, pool_size, flags, beta, 
                opt_type, total_score_saxs, total_score_cs, total_score_fret,
                total_score_jc, total_score_noe, total_score_pre, total_score_rdc,
                total_score_rh, bc_saxs, bc_cs, bc_fret, bc_alpha_vals_jc,
                bc_dist_vals_noe, bc_dist_vals_pre, bc_rdc, bc_rh
                )

        # calculate scores for unoptimized data types (JC, NOE, PRE)
        if not flags[pre_name]:
            opt_rmse_PRE, max_score_PREs = pre_optimization_ensemble(exp_data, bc_data, new_indices)[:2]
        if not flags[noe_name]:
            opt_rmse_NOE, max_score_NOEs = noe_optimization_ensemble(exp_data, bc_data, new_indices)[:2]
        if not flags[jc_name]:
            opt_rmse_JC, max_score_JC = jc_optmization_ensemble(exp_data, bc_data, new_indices)[:2]

        s = [
            i, accepted, max_score_SAXS, max_score_CS, max_score_FRET, 
            max_score_JC, max_score_NOEs, max_score_PREs, max_score_RDCs,
            max_score_RH, opt_rmse_SAXS, opt_rmse_CS, opt_rmse_FRET, 
            opt_rmse_JC, opt_rmse_NOE, opt_rmse_PRE, opt_rmse_RDC, 
            opt_rmse_RH
            ]
        final_results.append(s)
        final_indices.append(new_indices)
        final_best_jcoups.append(best_jcoups)
    
    return final_results, final_indices, final_best_jcoups

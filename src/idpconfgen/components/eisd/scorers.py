"""
Module contains all functions required to perform a single property
maximum loglike optimization.

Inspired/imported from:
* https://github.com/THGLab/X-EISD/blob/master/eisd/scorers.py
* https://github.com/Oufan75/X-EISD/blob/master/eisd/scorers.py
"""
import numpy as np

from idpconfgen.components.eisd import(
    saxs_name,
    cs_name,
    fret_name,
    jc_name,
    noe_name,
    pre_name,
    rdc_name,
    rh_name
    )


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


def calc_gamma(len_saxs, num_res, qmax, qmin):
    # Generalizable way to calculate gamma parameter
    # accournding to shannon theory
    N_q = len_saxs
    N_s = num_res * (qmax - qmin) / np.pi
    
    return N_s / N_q


def vect_calc_opt_params_jc(
    alpha1,
    alpha2,
    exp_j,
    exp_sig,
    mus,
    sigs,
    ):
    """
    Takes an average phi value and experimental data point
    calculates optimal A, B, C parameters
    """
    a = np.zeros((alpha1.shape[0], 3, 3))
    b = np.zeros((alpha1.shape[0], 3))

    a[:, 0, 0] = 1.0 / (sigs[0] ** 2.0) + ((alpha2 / exp_sig) ** 2.0)
    a[:, 1, 1] = 1.0 / (sigs[1] ** 2.0) + ((alpha1 / exp_sig) ** 2.0)
    a[:, 2, 2] = 1.0 / (sigs[2] ** 2.0) + 1.0 / (exp_sig ** 2.0)

    a[:, 0, 1] = alpha1 * alpha2 / (exp_sig ** 2.0)
    a[:, 1, 0] = alpha1 * alpha2 / (exp_sig ** 2.0)
    a[:, 0, 2] = alpha2 / (exp_sig ** 2.0)
    a[:, 2, 0] = alpha2 / (exp_sig ** 2.0)
    a[:, 1, 2] = alpha1 / (exp_sig ** 2.0)
    a[:, 2, 1] = alpha1 / (exp_sig ** 2.0)

    b[:, 0] = mus[0] / (sigs[0] ** 2.0) + exp_j * alpha2 / (exp_sig ** 2)
    b[:, 1] = mus[1] / (sigs[1] ** 2.0) + exp_j * alpha1 / (exp_sig ** 2)
    b[:, 2] = mus[2] / (sigs[2] ** 2.0) + exp_j / (exp_sig ** 2)

    opt_params = np.array([np.linalg.solve(a[i], b[i]) for i in range(a.shape[0])])  # shape: (47,3)

    return opt_params


def vect_calc_score_jc(
    alpha1,
    alpha2,
    exp_j,
    exp_sig,
    opt_params,
    mus,
    sigs
    ):
    """
    Calculates score for a single phi angle/residue, 
    can be used on a single protein or ensemble.
    
    Returns
    -------
    f : float
        Total score calculated
    
    f_comps : array
        [f_a, f_b, f_c, f_err]
    """
    f_a = normal_loglike(opt_params[:,0], mus[0], sigs[0])
    f_b = normal_loglike(opt_params[:,1], mus[1], sigs[1])
    f_c = normal_loglike(opt_params[:,2], mus[2], sigs[2])
    err = exp_j - opt_params[:,0] * alpha2 - opt_params[:,1] * alpha1 - opt_params[:,2]
    f_err = normal_loglike(err, 0 , exp_sig)
    f = f_a + f_b + f_c + f_err
    f_comps = [f_a, f_b, f_c, f_err]
    
    return f, f_comps

# ------------------------------------------------------------
# Below here we have the functions to score each ensemble
# by a specific experimental module
# ------------------------------------------------------------


def saxs_optimization_ensemble(
    exp_data,
    bc_data,
    indices,
    ens_size,
    nres,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    # prepare data
    exp_saxs = exp_data[saxs_name].data['value'].values
    exp_sigma = exp_data[saxs_name].data['error'].values
    
    if indices is None:
        bc_saxs = old_vals - \
            (bc_data[saxs_name].data.values[popped_structure, :] - \
            bc_data['saxs'].data.values[new_index, :] ) / ens_size
        
    else:
        bc_ensemble = bc_data[saxs_name].data.values[indices, :]
        bc_saxs = np.mean(bc_ensemble, axis=0)
        
    # optimization
    opt_params = calc_opt_params(bc_saxs, exp_saxs, exp_sigma, bc_data[saxs_name].sigma)
    
    g = calc_gamma(
        exp_data[saxs_name].data.shape[0],
        nres,
        np.max(exp_data[saxs_name].data['index']),
        np.min(exp_data[saxs_name].data['index'])
        )
    
    f, f_comps = calc_score(bc_saxs, exp_saxs, exp_sigma, bc_data[saxs_name].sigma, opt_params, gamma=g)

    error = (exp_saxs - bc_saxs) ** 2 
    rmse = np.mean(error) ** 0.5
    total_score = np.sum(f)
    
    return rmse, total_score, bc_saxs, error


def cs_optimization_ensemble(
    exp_data,
    bc_data,
    ens_size,
    indices,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    """
    Main logic for chemical shift scoring.

    Parameters
    ----------
    exp_data : dict
        Dictionary of experimental values of chemical shifts.
        First key-layer indicates type of CS (C, CA, CB, N, H, HA).
        For each key, there exists 2 values, the CS assignment [0]
        and the experimental error associated [1].
    
    bc_data : dict
        Dictionary of back-calculated values of chemical shifts.
        Format follows `exp_data`.
    
    indices : list
        List of indices of experimental and back-calculated data
        to include. Defaults to None.
        # !!! note here, may not be required based on new logic
    
    """
    # TODO: if incorrect shape, shave off values from bc_data
    # to match exp_data
    exp_cs = exp_data[cs_name].data['value'].values
    exp_sigma = exp_data[cs_name].data['error'].values
    atom_types = exp_data[cs_name].data['atomname'].values
    
    if indices is None:
        bc_cs = old_vals - \
            (bc_data[cs_name].data.values[popped_structure, :] - \
            bc_data['cs'].data.values[new_index, :]) / ens_size
    else:
        bc_ensemble = bc_data[cs_name].data.values[indices, :]
        bc_cs = np.mean(bc_ensemble, axis=0)
        
    bc_sigma = np.array([bc_data[cs_name].sigma[atom_type] for atom_type in atom_types])
    
    opt_params = calc_opt_params(bc_cs, exp_cs, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc_cs, exp_cs, exp_sigma, bc_sigma, opt_params)
    
    error = (exp_cs - bc_cs) ** 2.0
    rmse = np.mean(error) ** 0.5
    total_score = np.sum(f)
    
    return rmse, total_score, bc_cs, error


def fret_optimization_ensemble(
    exp_data,
    bc_data,
    ens_size,
    indices,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    """
    Main logic for fret scoring module
    """
    # prepare data
    exp = exp_data[fret_name].data
    exp_sigma = exp_data[fret_name].sigma

    if indices is None:
        bc = old_vals - \
            (bc_data[fret_name].data.values[popped_structure, :] - \
            bc_data['fret'].data.values[new_index, :]) / ens_size
    else:
        bc_ensemble = bc_data[fret_name].data.values[indices, :]
        bc = np.mean(bc_ensemble, axis=0)

    bc_sigma = bc_data[fret_name].sigma

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_sigma, opt_params)

    error = (exp - bc) ** 2.0
    rmse = np.mean(error) ** 0.5
    total_score = np.sum(f)

    return rmse, total_score, bc, error

    
def jc_optimization_ensemble(
    exp_data,
    bc_data,
    ens_size,
    indices,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    """
    Main logic for J-coupling scoring.
    """
    exp = exp_data[jc_name].data['value'].values
    exp_sigma = exp_data[jc_name].data['error'].values

    if indices is None:
        pop_alpha = bc_data[jc_name].data.values[popped_structure, :]
        add_alpha = bc_data[jc_name].data.values[new_index, :]

        bc_alpha1 = old_vals[0] - (pop_alpha - add_alpha) / ens_size
        bc_alpha2 = old_vals[1] - (np.square(pop_alpha) - np.square(add_alpha)) / ens_size
        
        assert bc_alpha1.shape == bc_alpha2.shape
    else:
        bc_ensemble_alpha1 = bc_data[jc_name].data.values[indices, :]
        bc_alpha1 = np.mean(bc_ensemble_alpha1, axis=0)

        bc_ensemble_alpha2 = np.square(bc_data[jc_name].data.values[indices, :])
        bc_alpha2 = np.mean(bc_ensemble_alpha2, axis=0)

    bc_sigma = [bc_data[jc_name].sigma[i] for i in ["A", "B", "C"]]
    bc_mu = [bc_data[jc_name].mu[i] for i in ["A", "B", "C"]]

    opt_params = vect_calc_opt_params_jc(
        bc_alpha1,
        bc_alpha2, 
        exp,
        exp_sigma,
        bc_mu,
        bc_sigma
        )
    
    f, f_comps = vect_calc_score_jc(
        bc_alpha1,
        bc_alpha2,
        exp,
        exp_sigma,
        opt_params,
        bc_mu,
        bc_sigma
        )
    
    error = (opt_params[:,0] * bc_alpha2 + opt_params[:,1] * bc_alpha1 + opt_params[:,2] - exp) ** 2.0
    rmse = np.mean(error) ** 0.5
    scores = np.sum(f)
    jcoup_vals = list(opt_params[:,0] * bc_alpha2 + opt_params[:,1] * bc_alpha1 + opt_params[:,2])

    return rmse, scores, jcoup_vals, [bc_alpha1, bc_alpha2]


def noe_optimization_ensemble(
    exp_data,
    bc_data,
    indices,
    ens_size,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    """
    Main logic for NOE scoring.
    """
    exp_distance = exp_data[noe_name].data['dist_value'].values
    upper_bound_value = exp_data[noe_name].data['upper'].values 
    lower_bound_value = exp_data[noe_name].data['lower'].values
    
    assert exp_distance.shape == upper_bound_value.shape == lower_bound_value.shape
    
    # uncertainty range should be made a parameter; confirm if compatible with exp file
    range_val = upper_bound_value + lower_bound_value
    exp_sigma = range_val / 2.0
    
    # load long range noe bc index
    # calculating inverse 6 average
    if indices is None:
        popped = np.power(bc_data[noe_name].data.values[popped_structure], -6.0)
        added = np.power(bc_data[noe_name].data.values[new_index], -6.0)
        # check that "ens_size" could be 100. ?
        avg_distance = (np.power(old_vals, -6.0) * ens_size - (popped - added) ) / ens_size
        avg_distance = np.power(avg_distance, (-1. / 6.))
    else:
        bc_ensemble = np.power(bc_data[noe_name].data.values[indices], -6.0)
        avg_distance = np.power(np.mean(bc_ensemble, axis=0), (-1. / 6.))
    
    # optimization
    opt_params = calc_opt_params(avg_distance, exp_distance, exp_sigma, bc_data[noe_name].sigma)
    f, f_comps = calc_score(avg_distance, exp_distance, exp_sigma, bc_data[noe_name].sigma, opt_params)
    
    error = (exp_distance - avg_distance) ** 2.0
    rmse = np.mean(error) ** 0.5
    total_score = np.sum(f)
    
    return rmse, total_score, avg_distance, error


def pre_optimization_ensemble(
    exp_data,
    bc_data,
    indices,
    ens_size,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    """
    Main logic for PRE scoring function.
    """
    # prepare data
    exp_distance = exp_data[pre_name].data['dist_value'].values
    upper_bound_value = exp_data[pre_name].data['upper'].values
    lower_bound_value = exp_data[pre_name].data['lower'].values
    
    assert exp_distance.shape == upper_bound_value.shape == lower_bound_value.shape
    
    range_val = upper_bound_value + lower_bound_value
    exp_sigma = range_val / 2.0
    
    if indices is None:
        popped = np.power(bc_data[pre_name].data.values[popped_structure, :], -6.0)
        added = np.power(bc_data[pre_name].data.values[new_index, :], -6.0)
        avg_distance = (np.power(old_vals, -6.0) * ens_size - (popped - added)) / ens_size
        avg_distance = np.power(avg_distance, (-1. / 6.))
    else:
        bc_ensemble = np.power(bc_data[pre_name].data.values[indices, :], -6.0)
        avg_distance = np.power(np.mean(bc_ensemble, axis=0), (-1./6.))
    
        # optimization
    opt_params = calc_opt_params(avg_distance, exp_distance, exp_sigma, bc_data[pre_name].sigma)
    f, f_comps = calc_score(avg_distance, exp_distance, exp_sigma, bc_data[pre_name].sigma, opt_params)
    
    error = (exp_distance - avg_distance) ** 2.0
    rmse = np.mean(error) ** 0.5
    total_score = np.sum(f)
    
    return rmse, total_score, avg_distance, error 


def rdc_optimization_ensemble(
    exp_data,
    bc_data,
    indices,
    ens_size,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    # prepare data
    exp = exp_data[rdc_name].data['value'].values
    exp_sigma = exp_data[rdc_name].data['error'].values
    
    assert exp.shape == exp_sigma.shape

    if indices is None:
        bc = old_vals - \
            (bc_data[rdc_name].data.values[popped_structure, :] - \
            bc_data[rdc_name].data.values[new_index, :]) / ens_size
    else:
        bc_ensemble = bc_data[rdc_name].data.values[indices, :]
        bc = np.mean(bc_ensemble, axis=0)

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_data[rdc_name].sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_data[rdc_name].sigma, opt_params)
    
    error = (exp - bc) ** 2.0
    rmse = np.mean(error) ** 0.5
    total_score = np.sum(f)

    return rmse, total_score, bc, error


def rh_optimization_ensemble(
    exp_data,
    bc_data,
    indices,
    ens_size,
    old_vals=None,
    popped_structure=None,
    new_index=None,
    ):
    # prepare data
    exp = exp_data[rh_name].data
    exp_sigma = exp_data[rh_name].sigma

    if indices is None:
        bc = old_vals - \
            (bc_data[rh_name].data.values[popped_structure, :] - \
            bc_data[rh_name].data.values[new_index, :]) / ens_size
    else:
        bc_ensemble = bc_data[rh_name].data.values[indices, :]
        bc = np.mean(bc_ensemble, axis=0)

    bc_sigma = bc_data[rh_name].sigma

    # optimization
    opt_params = calc_opt_params(bc, exp, exp_sigma, bc_sigma)
    f, f_comps = calc_score(bc, exp, exp_sigma, bc_sigma, opt_params)
    
    error = (exp - bc) ** 2.0
    rmse = np.mean(error)**0.5
    total_score = np.sum(f)

    return rmse, total_score, bc[0], error

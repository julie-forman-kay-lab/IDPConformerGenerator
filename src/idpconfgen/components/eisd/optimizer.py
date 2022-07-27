"""
Inspired and imported from 
https://github.com/THGLab/X-EISD/blob/master/eisd/optimizer.py
https://github.com/Oufan75/X-EISD/blob/master/eisd/optimizer.py
"""
import numpy as np

from idpconfgen.components.eisd import (
    modes,
    eisd_run_all,
    opt_max,
    opt_mc,
    )
from idpconfgen.components.eisd.scorers import *


def monte_carlo(beta, old_total_score, new_total_score):
    new_probability = np.exp(beta * new_total_score)
    old_probability = np.exp(beta * old_total_score)
    
    # to deal with runtime error caused by large score values
    if np.any(np.isinf([old_probability, new_probability])):
        # print('Runtime error... reset beta value')
        beta = 500. / new_probability
        new_probability = np.exp(beta * new_total_score)
        old_probability = np.exp(beta * old_total_score)
        
    # accept criterion
    return np.random.random_sample() < min(1, new_probability/old_probability)


class XEISD(object):
    """
    This is the API to the X-EISD scoring to calculate and/or optimize log-likelihood of a
    disordered protein ensemble.
    
    Parameters
    ----------
    exp_data : dict
        Experimental data files with uncertainties.
    bc_data : dict
        Back calculation files with uncertainties.
    ens_size : int
        number of candidate conformers.
    """
    def __init__(self, exp_data, bc_data, nres, ens_size=None):
        self.exp_data = exp_data
        self.bc_data = bc_data
        self.resnum = nres
        self.ens_size = ens_size
        if ens_size is None:
            # If pool size not provided. Use the JC back-calc size
            self.ens_size = bc_data[jc_name].data.shape[0]


    def calc_scores(self, dtypes, ens_size, indices=None):
        '''
        Parameters
        ----------
        dtypes : list
            list of data types to score
        
        ens_size : int
            Only used when indices not specified.
            Used to randomly select subset to score.
        
        indices : ndarray, optional
            This is the fastest way to get the EISD score and RMSD of selected properties for a given set of indices.
            shape: (size_of_ensemble, )
            Defaults to None.
        
        Return
        ----------
        scores : dict 
            Dictionary of RMSDs.
            X-EISD scores and ensemble averages for selected conformers.
        '''
        if indices is None:
            indices = np.random.choice(np.arange(self.ens_size), ens_size, replace=False)

        # initiate dict to store scores
        scores = {}
            
        for name in self.exp_data.keys():
            scores[name] = [0, 0, 0]
            if name == jc_name:
                scores[name] = [0, 0, 0, [0]]
        
        # The first 3 items list[:3] is taken as the "rmse", "score", "back-calculation"
        # last "error" return is not used for X-EISD
        if jc_name in dtypes:
            scores[jc_name] = list(jc_optimization_ensemble(self.exp_data, self.bc_data, indices))
        if saxs_name in dtypes:
            scores[saxs_name] = list(saxs_optimization_ensemble(self.exp_data, self.bc_data, indices, nres=self.resnum))[:3]
        if cs_name in dtypes:
            scores[cs_name] = list(cs_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if fret_name in dtypes:
            scores[fret_name] = list(fret_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if noe_name in dtypes:
            scores[noe_name] = list(noe_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if pre_name in dtypes:
            scores[pre_name] = list(pre_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if rdc_name in dtypes:
            scores[rdc_name] = list(rdc_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
        if rh_name in dtypes:
            scores[rh_name] = list(rh_optimization_ensemble(self.exp_data, self.bc_data, indices))[:3]
     
        return scores


    def optimize(
        self,
        epochs,
        ens_size,
        opt_type=opt_max,
        mode=eisd_run_all,
        beta=0.1,
        iters=100,
        ):
        """
        Parameters
        ----------        
        opt_type : str
            The optimization type should be 'mc' or 'max', 'mc' for Metropolis Monte Carlo, 
            Defaults to 'max' for score maximization method.
        
        epochs : int
            Number of optimization trials.
        
        ens_size : int
            Number of conformers in the ensemble.
            
        mode : str or list
            Data types to optimize.
            Defaults to "all".
            
        beta : float
            Temperature parameter for MC optimization.
            Defaults to 0.1.
            
        iters : int
            Number of conformer exchange attempts.
            Defaults to 100.
            
        Returns
        -------
        final_results : list
        
        final_indices : list
        
        final_best_jcoups : list
        
        """
        # switch the property
        flags = modes(mode, self.exp_data.keys())

        final_results = []
        final_indices = []
        final_best_jcoups = []

        for it in range(epochs):
            # initial scores
            indices = list(np.random.choice(np.arange(self.ens_size), ens_size, replace=False))
            old_scores = self.calc_scores([key for key in flags if flags[key]], indices)    

            new_scores = {}
            for name in flags:
                new_scores[name] = [0, 0, 0]
                if name == jc_name:
                    new_scores[name] = [0, 0, 0, [0]]
            accepted = 0
            
            for iterations in range(iters):
                pop_index = np.random.randint(0, ens_size, 1)[0]
                popped_structure = indices[pop_index]
                indices.pop(pop_index)
                struct_found = False
                while not struct_found:
                    new_index = np.random.randint(0, self.ens_size, 1)[0]
                    if new_index != popped_structure and new_index not in indices:
                        indices.append(new_index)
                        struct_found = True

                for prop in flags:
                    if flags[prop]:
                        if prop == saxs_name:
                            new_scores[saxs_name] = \
                                list(saxs_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[saxs_name][2], popped_structure,
                                new_index, ens_size, self.resnum))[:3]

                        if prop == cs_name:
                            new_scores[cs_name] = \
                                list(cs_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[cs_name][2], popped_structure,
                                new_index, ens_size))[:3]

                        if prop == fret_name:
                            new_scores[fret_name] = \
                                list(fret_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[fret_name][2], popped_structure,
                                new_index, ens_size))[:3]

                        if prop == jc_name:
                            new_scores[jc_name] = \
                                list(jc_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[jc_name][3], popped_structure,
                                new_index, ens_size))

                        if prop == noe_name:
                            new_scores[noe_name] = \
                                list(noe_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[noe_name][2], popped_structure,
                                new_index, ens_size))[:3]
 
                        if prop == pre_name:
                            new_scores[pre_name] = \
                                list(pre_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[pre_name][2], popped_structure,
                                new_index, ens_size))[:3]

                        if prop == rdc_name:
                            new_scores[rdc_name] = \
                                list(rdc_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[rdc_name][2], popped_structure,
                                new_index, ens_size))[:3]
                                
                        if prop == rh_name:
                            new_scores[rh_name] \
                                = list(rh_optimization_ensemble(self.exp_data, 
                                self.bc_data, None, old_scores[rh_name][2], popped_structure,
                                new_index, ens_size))[:3]
            
                old_total_score = np.sum([old_scores[key][1] for key in old_scores])
                new_total_score = np.sum([new_scores[key][1] for key in new_scores])

                # optimization
                if opt_type == opt_max:
                    to_accept = old_total_score < new_total_score
                elif opt_type == opt_mc:
                    to_accept = monte_carlo(beta, old_total_score, new_total_score)
            
                if not to_accept:
                    indices.pop(-1)
                    indices.append(popped_structure)
                else:
                    for prop in flags:
                        old_scores[prop] = new_scores[prop]

                    accepted = accepted + 1         

            s = [it, accepted]
            for prop in flags:
                # calculate scores for unoptimized data types
                if not flags[prop]:
                    if prop == pre_name:
                        old_scores[pre_name][:2] = \
                            pre_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == jc_name:
                        old_scores[jc_name][:2] = \
                            jc_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == cs_name:
                        old_scores[cs_name][:2] = \
                            cs_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == fret_name:
                        old_scores[fret_name][:2] = \
                            fret_optimization_ensemble(self.exp_data, self.bc_data, indices)[:2]
                    if prop == saxs_name:
                        old_scores[saxs_name][:2] = \
                            saxs_optimization_ensemble(self.exp_data, self.bc_data, indices, nres=self.resnum)[:2]  # noqa: E501
                # aggregate results
                s.extend(old_scores[prop][:2])

            final_results.append(s)
            final_indices.append(indices)
            final_best_jcoups.append(old_scores[jc_name][2])
        
        return final_results, final_indices, final_best_jcoups

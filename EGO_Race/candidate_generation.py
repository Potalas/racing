from BayesOpt import BO
from BayesOpt.BayesOpt import generate_candidate_global
from BayesOpt.SearchSpace import ContinuousSpace, OrdinalSpace, ProductSpace, NominalSpace
from BayesOpt.Surrogate import RandomForest
from BayesOpt.base import Solution

from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)





class modified_BO(BO):
    def __init__(self, params, try_use_gp = False, **kwargs):
        #Extract different parameter types
        param_names = [x for (x,y) in zip(params['names'], params['isFixed']) if not y]
        cont_params = [x for (x,y) in zip(param_names, params['types']) if y == "r"]
        ordinal_params = [x for (x,y) in zip(param_names, params['types']) if y == "i"]
        nominal_params = [x for (x,y) in zip(param_names, params['types']) if y == "c" or y == "o"]
        fixed_params = [x for (x,y) in zip(params['names'], params['isFixed']) if y]
        self.fixed_params = {x : params['domain'][x] for x in fixed_params}
        self.param_names = cont_params + ordinal_params + nominal_params
        
        #Turn into searchspace
        search_space = None
        if len(cont_params) > 0:
            search_space = ContinuousSpace([params['domain'][x] for x in cont_params], name=cont_params)
        if len(ordinal_params) > 0:
            search_space_ordinal = OrdinalSpace([params['domain'][x] for x in ordinal_params], name=ordinal_params)
            if search_space is None:
                search_space = search_space_ordinal
            else:
                search_space += search_space_ordinal
        if len(nominal_params) > 0:
            search_space_nominal = NominalSpace([params['domain'][x] for x in nominal_params], name=nominal_params)    
            if search_space is None:
                search_space = search_space_nominal
            else:
                search_space += search_space_nominal
        if try_use_gp and len(ordinal_params) + len(nominal_params) == 0:
            from GaussianProcess.gpr import GaussianProcess
            from GaussianProcess.trend import constant_trend
            
            dim = len(cont_params)
            mean = constant_trend(dim, beta=0)    
            lb = np.array([x[0] for x in search_space.bounds])
            ub = np.array([x[1] for x in search_space.bounds])
            print(lb,ub)
            thetaL = 1e-10 * (ub - lb) * np.ones(dim)
            thetaU = 2 * (ub - lb) * np.ones(dim)
            theta0 = np.random.rand(dim) * (thetaU - thetaL) + thetaL
            print(thetaL, thetaU)
            model = GaussianProcess(mean=mean, corr='matern',
                        theta0=theta0, thetaL=thetaL, thetaU=thetaU,
                        nugget=1e-10, noise_estim=False,
                        optimizer='BFGS', wait_iter=5, random_start=30 * dim,
                        likelihood='concentrated', eval_budget=150 * dim)
            super(modified_BO, self).__init__(search_space, None, model, max_eval = 100, infill = "MGFI", **kwargs, t0 = 2)

        else:
#             super(modified_BO, self).__init__(search_space, None, RandomForest(max_features = 'sqrt', min_samples_leaf = 1, n_estimators = 1000), 
#                                               max_eval = 100, infill = "MGFI", **kwargs, t0 = 2)
            levels = search_space.levels if hasattr(search_space, 'levels') else None
            super(modified_BO, self).__init__(search_space, None, RandomForest(levels = levels), max_eval = 100, infill = "MGFI", **kwargs, t0 = 2)
        

    def add_data(self, dt):
        DOE = Solution(np.array(dt[self.param_names]), var_name = self.param_names, n_obj=1, fitness=np.array(dt['val']), n_eval = 1)
        DOE = self.after_eval_check(DOE)

        if hasattr(self, 'data'):
            #If new dt contains all previous data as subset, overwrite, otherwise add
            if np.array([DOE[:len(self.data)] == self.data]).all():
                self.data = DOE
            else:
                self.data += DOE
        else:
            self.data = DOE
        self.fit_and_assess()
        self.plugin = min(self.surrogate.predict(self.data))
        
    def add_surrogate_data(self, dt):
        DOE = Solution(np.array(dt[self.param_names]), var_name = self.param_names, n_obj=1, n_eval = 1)
        def pred_func(s):
            s = self.frange * self.surrogate.predict(s.reshape(1,-1))[0] + self.fmin
            return s
        fs = list(map(pred_func, DOE))
        DOE.fitness = fs
        DOE = self.after_eval_check(DOE)
        self.data += DOE
        self.fit_and_assess()
        self.plugin = min(self.surrogate.predict(self.data))
    

    def get_candidates(self, nr_candidates, seed = None):
        ignore_acquisition = False
        if seed is not None:
            if seed < 0:
                ignore_acquisition = True
                np.random.seed(int(-1*seed))
            else:
                np.random.seed(seed)
#         print("get_candidates fct")
        self.n_point = nr_candidates
        if not hasattr(self, 'data'):
            candidates = self._space.sampling(nr_candidates)
        else:
            candidates, values = self.arg_max_acquisition2(ignore_acquisition)
#         print("done with acq")
        df = pd.DataFrame.from_records(candidates, columns=self.param_names)
#         print(values)
        for k,v in self.fixed_params.items():
            df[k] = v
#         print("returning")
        return df
    

    def arg_max_acquisition2(self, ignore_acquisition = False):
        """
        Global Optimization of the acqusition function / Infill criterion
        Returns
        -------
            candidates: tuple of list,
                candidate solution (in list)
            values: tuple,
                criterion value of the candidate solution
        """
        dx = True if self._optimizer == 'BFGS' else False
        
        if self.n_job > 1 and self.n_point > 1:
            # TODO: fix this issue once for all!
#             try:  
#                 self.pool.restart() # restart the pool in case it is terminated before
#             except AssertionError:
#                 pass
            t_values = [np.exp(self.t * np.random.randn()) for i in range(self.n_point)]
            print(t_values)
            generate_candidate_partial = partial(generate_candidate_global, plugin = self.plugin, surrogate = self.surrogate, 
                                                 minimize = self.minimize, _random_start = self._random_start, 
                                                 _space = self._space, eq_func = self.eq_func, ineq_func = self.ineq_func,
                                                 _eval_type = self._eval_type, _wait_iter = self._wait_iter, _max_eval = self._max_eval)
#             __ = self.pool.map(generate_candidate_partial, t_values)
            p = Pool(min([self.n_job, len(t_values), cpu_count()]))
            __ = p.map(generate_candidate_partial, t_values)
            p.close()
        else:
            t_value = np.exp(self.t * np.random.randn())
            __ = [generate_candidate_global(t_value, plugin = self.plugin, surrogate = self.surrogate, 
                                            minimize = self.minimize, _random_start = self._random_start, 
                                            _space = self._space, eq_func = self.eq_func, ineq_func = self.ineq_func,
                                            _eval_type = self._eval_type, _wait_iter = self._wait_iter, _max_eval = self._max_eval,
                                           ignore_acquisition = ignore_acquisition)]
#         else:
#             criteria = [self._acquisition(plugin, dx=dx) for i in range(self.n_point)]
#             __ = [list(self._argmax_multistart(_)) for _ in criteria]

        candidates, values = tuple(zip(*__))

        return candidates, values
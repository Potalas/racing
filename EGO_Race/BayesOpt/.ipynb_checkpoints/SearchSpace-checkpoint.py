"""
@author: Hao Wang
@email: wangronin@gmail.com
"""
from pdb import set_trace

import six, json
from copy import copy, deepcopy

import numpy as np
from numpy.random import randint, rand
from abc import abstractmethod
from pyDOE import lhs

class SearchSpace(object):
    def __init__(self, bounds, var_name, name):
        """Search Space Base Class

        Parameters
        ---------
        bounds : (list of) list,
            lower and upper bound for continuous/ordinal parameter type
            categorical values for nominal parameter type.
            The dimension of the space is determined by the length of the 
            nested list
        var_name : (list of) str,
            variable name per dimension. If only a string is given for multiple 
            dimensions, variable names are created by appending counting numbers
            to the input string. 
        name : str,
            search space name. It is typically used as the grouping variable
            when converting the Solution object to dictionary, allowing for 
            vector-valued search parameters. See 'to_dict' method below.

        Attributes
        ----------
        dim : int,
            dimensinality of the search space
        bounds : a list of lists,
            each sub-list stores the lower and upper bound for continuous/ordinal variable
            and categorical values for nominal variable
        levels : a list of lists,
            each sub-list stores the categorical levels for every nominal variable. It takes
            `None` value when there is no nomimal variable
        precision : a list of double,
            the numerical precison (granularity) of continuous parameters, which usually  
            very practical in real-world applications
        var_name : a list of str,
            variable names per dimension 
        var_type : a list of str, 
            variable type per dimension, 'C': continuous, 'N': nominal, 'O': ordinal
        C_mask : a bool array,
            the mask array for continuous variables
        O_mask : a bool array,
            the mask array for integer variables
        N_mask : a bool array,
            the mask array for discrete variables
        id_C : an int array,
            the index array for continuous variables
        id_O : an int array,
            the index array for integer variables
        id_N : an int array,
            the index array for discrete variables
        """
        if hasattr(bounds[0], '__iter__') and not isinstance(bounds[0], str):
            self.bounds = [tuple(b) for b in bounds]
        else:
            self.bounds = [tuple(bounds)]
            
        self.dim = len(self.bounds)
        self.name = name
        self.var_type = None 
        self.levels = None
        self.precision = None

        if var_name is not None:
            if isinstance(var_name, str):
                if self.dim > 1:
                    var_name = [var_name + '_' + str(_) for _ in range(self.dim)]
                else:
                    var_name = [var_name]
            assert len(var_name) == self.dim
            self.var_name = var_name

    @abstractmethod
    def sampling(self, N=1):
        """The output is a list of shape (N, self.dim)
        """
        pass
    
    def _set_index(self):
        self.C_mask = np.asarray(self.var_type) == 'C'  # Continuous
        self.O_mask = np.asarray(self.var_type) == 'O'  # Ordinal
        self.N_mask = np.asarray(self.var_type) == 'N'  # Nominal 
        
        self.id_C = np.nonzero(self.C_mask)[0]
        self.id_O = np.nonzero(self.O_mask)[0]
        self.id_N = np.nonzero(self.N_mask)[0]

    def _set_levels(self):
        """Set categorical levels for all nominal variables
        """
        if hasattr(self, 'id_N') and len(self.id_N) > 0:
            self.levels = {i : self.bounds[i] for i in self.id_N}
            self._n_levels = {i : len(self.bounds[i]) for i in self.id_N}
        else:
            self.levels, self._n_levels = None, None 

    # TODO: this functionality should be more generic. See `pandas.to_dict`
    def to_dict(self, solution):
        """Convert a `Solution` object to a dictionary

        Parameters
        ----------
        solution : Solution
            an instance of .base.Solution class

        Returns
        -------
        dict
            A dictionary containing the data of solution
        """
        if self.name is None:
            return solution.to_dict() 
        else:
            return dict(zip(self.name, solution.tolist()))

    def round(self, X):
        if self.precision is not None:
            X = deepcopy(X)
            if not isinstance(X[0], list):
                X = [X]
            for k, v in self.precision.items():
                for i in range(len(X)):
                    X[i][k] = np.round(X[i][k], v)
        return X

    def __len__(self):
        return self.dim

    def __iter__(self):
        pass

    def __add__(self, space):
        """Direct Sum of two Spaces
        """
        assert isinstance(space, SearchSpace)
        return ProductSpace(self, space)

    def __radd__(self, space):
        return self.__add__(space)

    def __mul__(self, N):
        """Replicate a Space N times
        """
        N = int(N)
        s = deepcopy(self)
        s.dim = int(self.dim * N)
        s.var_type *= N
        s.bounds *= N
        s.var_name = ['{}_{}'.format(v, k) for k in range(N) for v in self.var_name]
        return s

    def __rmul__(self, N):
        return self.__mul__(N)
    

# TODO: maybe implement the scalar multiplication for ProductSpace
class ProductSpace(SearchSpace):
    """Cartesian product of the search spaces
    """
    def __init__(self, spaceL, spaceR):
        # setup the space names
        nameL = spaceL.name if isinstance(spaceL, ProductSpace) else [spaceL.name] 
        nameR = spaceR.name if isinstance(spaceR, ProductSpace) else [spaceR.name]
        self.name = nameL + nameR
        self.dim = spaceL.dim + spaceR.dim
        # TODO: check coincides of variable names
        self.var_name = spaceL.var_name + spaceR.var_name             
        self.bounds = spaceL.bounds + spaceR.bounds
        self.var_type = spaceL.var_type + spaceR.var_type

        self._subspaceL = deepcopy(spaceL)
        self._subspaceR = deepcopy(spaceR)
        self._set_index()
        self._set_levels()

        precision = copy(spaceL.precision) if spaceL.precision is not None else {}
        _ = {(k + spaceL.dim) : v for k, v in spaceR.precision.items()} \
            if spaceR.precision is not None else {}
        precision.update(_)
        self.precision = precision if precision else None 
    
    def sampling(self, N=1, method='uniform'):
        a = self._subspaceL.sampling(N, method)
        b = self._subspaceR.sampling(N, method)
        return [a[i] + b[i] for i in range(N)]
    
    def to_dict(self, solution):
        """Save a Solution instance to a dictionary 
        
        The result is grouped by sub-spaces, which is meant for vector-valued 
        parameters for the configuration 

        Parameters
        ----------
        solution : .base.Solution
            A solution object

        Returns
        -------
        dict
        """
        id1 = list(range(self._subspaceL.dim))
        id2 = list(range(self._subspaceL.dim, self.dim))
        L = solution[id1] if len(solution.shape) == 1 else solution[:, id1]
        R = solution[id2] if len(solution.shape) == 1 else solution[:, id2]
        return {**self._subspaceL.to_dict(L), **self._subspaceR.to_dict(R)}

    def __mul__(self, space):
        raise ValueError('Unsupported operation')

    def __rmul__(self, space):
        raise ValueError('Unsupported operation')


class ContinuousSpace(SearchSpace):
    """Continuous (real-valued) Search Space
    """
    def __init__(self, bounds, var_name='r', name=None, precision=None):
        super(ContinuousSpace, self).__init__(bounds, var_name, name)
        self.var_type = ['C'] * self.dim
        self._bounds = np.atleast_2d(self.bounds).T
        self._set_index()
        
        # set up precisions for each dimension
        if hasattr(precision, '__iter__'):
            assert len(precision) == self.dim
            self.precision = {i : precision[i] \
                for i in range(self.dim) if precision[i] is not None}
        else:
            if precision is not None:
                self.precision = {i : precision for i in range(self.dim)}
        
        assert all(self._bounds[0, :] < self._bounds[1, :])

    def __mul__(self, N):
        s = super(ContinuousSpace, self).__mul__(N)
        s._bounds = np.tile(s._bounds, (1, N))
        s._set_index()
        if self.precision is None:
            s.precision = None
        else:
            s.precision = {}
            for i in range(N):
                s.precision.update(
                    {(k + self.dim * i) : v for k, v in self.precision.items()}
                )
        return s
    
    def sampling(self, N=1, method='uniform'):
        lb, ub = self._bounds
        if method == 'uniform':   # uniform random samples
            X = ((ub - lb) * rand(N, self.dim) + lb).tolist()
        elif method == 'LHS':     # Latin hypercube sampling
            if N == 1:
                X = ((ub - lb) * rand(N, self.dim) + lb).tolist()
            else:
                X = ((ub - lb) * lhs(self.dim, samples=N, criterion='cm') + lb).tolist()

        return self.round(X)


class NominalSpace(SearchSpace):
    """Nominal (discrete) Search Space
    """
    def __init__(self, levels, var_name='d', name=None):
        levels = self._get_unique_levels(levels)
        super(NominalSpace, self).__init__(levels, var_name, name)
        self.var_type = ['N'] * self.dim
        self._levels = [np.array(b, dtype=object) for b in self.bounds]
        self._set_index()
        self._set_levels()
    
    def _get_unique_levels(self, levels):
        index = list(hasattr(l, '__iter__') and not isinstance(l, str) for l in levels)
        if any(index):
            return [
                list(set(levels[k] if i else [levels[k]])) \
                    for k, i in enumerate(index)
            ]
        else: 
            return [list(set(levels))]

    def __mul__(self, N):
        s = super(NominalSpace, self).__mul__(N)
        s._set_index()
        s._set_levels()
        return s
    
    def sampling(self, N=1, method='uniform'):
        # NOTE: `LHS` sampling does not apply here since nominal variable is not ordered
        res = np.empty((N, self.dim), dtype=object)
        for i in range(self.dim):
            idx = randint(0, self._n_levels[i], N)
            res[:, i] = [self.levels[i][_] for _ in idx]

        return res.tolist()


class OrdinalSpace(SearchSpace):
    """Ordinal (integer) Search Space
    """
    def __init__(self, bounds, var_name='i', name=None):
        super(OrdinalSpace, self).__init__(bounds, var_name, name)
        self.var_type = ['O'] * self.dim
        self._lb, self._ub = zip(*self.bounds)        # for sampling
        assert all(np.array(self._lb) < np.array(self._ub))
        self._set_index()

    def __mul__(self, N):
        s = super(OrdinalSpace, self).__mul__(N)
        s._lb, s._ub = s._lb * N, s._ub * N
        s._set_index()
        return s
    
    def sampling(self, N=1, method='uniform'):
        # TODO: adding LHS sampling here
        res = np.zeros((N, self.dim), dtype=int)
        for i in range(self.dim):
            res[:, i] = list(map(int, randint(self._lb[i], self._ub[i], N)))
        return res.tolist()


def from_dict(param, space_name=True):
    """Create a search space object from input dictionary

    Parameters
    ----------
    param : dict
        A dictionary that describes the search space
    space_name : bool, optional
        Whether a (multi-dimensional) subspace should be named. If this named space 
        is a subspace a whole search space, for a solution sampled from the whole space, its 
        components pertaining to this subspace will be grouped together under the key 
        `space_name`, when this solution is converted to a dictionary/json
        (see `SearchSpace.to_dict`).

    Returns
    -------
    SearchSpace
    """
    assert isinstance(param, dict)
    # construct the search space
    for i, (k, v) in enumerate(param.items()):
        bounds = v['range']
        if not hasattr(bounds[0], '__iter__') or isinstance(bounds[0], str):
            bounds = [bounds]

        bounds *= v['N']
        name = k if space_name else None

        # IMPORTANT: name argument is necessary for the variable grouping
        if v['type'] == 'r':        # real-valued parameter
            precision = v['precision'] if 'precision' in v else None 
            space_ = ContinuousSpace(
                bounds, var_name=k, name=name, 
                precision=precision
            ) 
        elif v['type'] == 'i':      # integer-valued parameter
            space_ = OrdinalSpace(bounds, var_name=k, name=name)
        elif v['type'] == 'c':      # category-valued parameter
            space_ = NominalSpace(bounds, var_name=k, name=name) 
        
        if i == 0:
            space = space_
        else:
            space += space_
    
    return space

def from_json(file):
    """Create a seach space from a json file

    Parameters
    ----------
    file : str
        Path to the input json file

    Returns
    -------
    SearchSpace
        an `SearchSpace` object converted from the json file
    """
    with open(file, 'r') as f:
        space = from_dict(json.load(f))
    return space

"""Create a seach space from a dict originating from irace

    This function expects a dictionary, which can be created from
    the irace parameters object by using the reticulate library
    This functions does not deal with conditional parameters yet,
    so these should be removed before calling this.

    Parameters
    ----------
    params : dict
        The dictionary containing the parameters

    Returns
    -------
    SearchSpace
        an `SearchSpace` object converted from the dictionary
    """
def from_irace(params):
    param_names = params['names']
    cont_params = [x for (x,y) in zip(param_names, params['types']) if y == "r"]
    ordinal_params = [x for (x,y) in zip(param_names, params['types']) if y == "i"]
    nominal_params = [x for (x,y) in zip(param_names, params['types']) if y == "c" or y == "o"]
#     fixed_params = [x for (x,y) in zip(params['names'], params['isFixed']) if y]
#     self.fixed_params = {x : params['domain'][x] for x in fixed_params}
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
    return search_space

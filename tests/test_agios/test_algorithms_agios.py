import numpy as np
import pandas as pd
import pytest

import corna.helpers as hl
from corna.algorithms import matrix_calc as algo


iso_tracer = ['C13']

no_atom_tracer = 2

eleme_corr = {}

na_dict = {'C': [0.95, 0.05], 'H':[0.98,0.01,0.01], 'S': [0.922297, 0.046832, 0.030872], 'O':[0.95,0.03,0.02], 'N': [0.8, 0.2]}

df = pd.DataFrame({'Name': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'}, \
   'Parent': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'}, \
    'Label': {0: 'C13_0', 1: 'C13_1', 2: 'C13_2'}, \
     'Intensity': {0: 0.3624, 1: 0.040349999999999997, 2: 0.59724999999999995}, \
      'Formula': {0: 'H4C2O2', 1: 'H4C2O2', 2: 'H4C2O2'}, \
      'info2': {0: 'culture_1', 1: 'culture_1', 2: 'culture_1'}, \
       'Sample': {0: 'sample_1', 1: 'sample_1', 2: 'sample_1'}})

correc_inten_dict = {'sample_1': {(0, 1): np.array([ 0.0619]), (0, 0): np.array([ 0.2456]), \
     (3, 0): np.array([ 0.0003]), (3, 1): np.array([ 0.0001]), (2, 1): np.array([ 0.0015]), \
      (2, 0): np.array([ 0.0071]), (5, 0): np.array([ 0.]), (5, 1): np.array([ 0.60045]), \
      (1, 0): np.array([ 0.06665]), (4, 1): np.array([ 0.]), (1, 1): np.array([ 0.0164]), \
       (4, 0): np.array([ 0.])}}

label_list = [(0, 1), (0, 0), (3, 0), (3, 1), (2, 1), (2, 0), (5, 0), \
                 (5, 1), (1, 0), (4, 1), (1, 1), (4, 0)]


def test_make_expected_na_matrix():
    assert algo.make_expected_na_matrix(0, [0, 1]) == np.array([1.])


def test_corr_matrix():
    iso_tracer = 'C'
    with pytest.raises(KeyError):
        c_matrix = algo.make_correction_matrix(iso_tracer, {'H': 1}, na_dict, ['H'])

def test():
    na = {'C': [0.95, 0.05],'H':[0.98,0.01,0.01], 'S': [0.922297, 0.046832, 0.030872], 'O17':[0.95,0.03], 'O18':[0.95, 0, 0.02], 'N': [0.8, 0.2]}
    iso_tracer = 'C'
    c_matrix = algo.make_correction_matrix(iso_tracer, {'C': 2,'H': 1, 'O':3}, na_dict, ['O'])
    print c_matrix

def test_make_all_corr_matrices():
    iso_tracer = 'C'
    with pytest.raises(KeyError):
        algo.make_all_corr_matrices(iso_tracer, {'H': 1}, na_dict, {'C': []})

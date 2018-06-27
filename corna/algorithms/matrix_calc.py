
import numpy as np
from numpy.linalg import pinv
import pandas as pd

from corna.constants import ISOTOPE_NA_MASS, KEY_ELE
from corna.helpers import get_isotope_element

def make_expected_na_matrix(N, pvec):
    """for a single labeled element, create the matrix M
    such that Mx=y where x is the actual distribution of input labels
    and y is the expected distribution of intensities with natural abundance

    N: number of atoms of this element
    pvec: expected isotopic distribution (e.g. [0.99,0.01])"""

    max_label=1+(N*(len(pvec)-1))
    correction_matrix = np.zeros((max_label,N+1))
    for i in range(N+1):
        column = np.zeros(i+1)
        column[-1]=1.0
        for nb in range(N-i):
            column = np.convolve(column, pvec)
        column.resize(max_label)
        correction_matrix[:, i] = column
    return correction_matrix

def add_indistinguishable_element(M,n,pvec):
    """to a matrix M formed by make_expected_na_matrix, add additional expected
    intensity corresponding to natural abundance from elements with same isotopic mass shift

    M: previous matrix formed by make_expected_na_matrix
    n: number of atoms of new element
    pvec: expected isotopic distribution of new element (e.g. [0.99,0.01])"""
    max_label = (n*(len(pvec)-1))
    M_new = np.zeros((M.shape[0]+max_label, M.shape[1]))
    M_new[:M.shape[0], :] = M
    for i in range(M.shape[1]):
        for j in range(n):
            M_new[:, i] = np.convolve(M_new[:, i], pvec)[:M_new.shape[0]]
    #print M_new
    return M_new

def make_correction_matrix(trac_atom, formuladict, na_dict, indist_elems):
    """create matrix M such that Mx=y where y is the observed isotopic distribution
    and x is the expected distribution of input labels

    atom_bag: dict of element:number of atoms in molecule (e.g. {'C':2,'O':1,'H':6})
    label_elem: element with input labeling
    indist_elems: elements with identical mass shift
    na_dict: dict of element:expected isotopic distribution
    :TODO This function relates to issue NCT-247. Need to change the function
    in more appropriate way.
    """
    M = make_expected_na_matrix(formuladict.get(trac_atom, 0), na_dict[trac_atom])
    for e in indist_elems:
        if e in formuladict:
            e1 = ISOTOPE_NA_MASS[KEY_ELE][e]
            M = add_indistinguishable_element(M, formuladict[e1], na_dict[e])
    return pinv(M)


def make_all_corr_matrices(isotracers, formula_dict, na_dict, eleme_corr):
    """
    This function forms correction matrix according to each isotracer.
    """
    corr_mats = {}
    for isotracer in isotracers:
        trac_atom = get_isotope_element(isotracer)
        try:
            indist_list = eleme_corr[trac_atom]
        except KeyError:
            indist_list = []
        corr_mats[isotracer] = make_correction_matrix(trac_atom, formula_dict, na_dict, indist_list)
    return corr_mats


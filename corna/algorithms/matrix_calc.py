
from copy import copy

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
    return M_new

def add_indistinguishable_element_for_autodetect(M, n,pvec, corr_limit):
    max_label = (n*(len(pvec)-1))
    M_new = np.zeros((M.shape[0], M.shape[1]))
    M_new[:M.shape[0], :] = np.eye(M.shape[1])
    for i in range(M.shape[1]):
        for j in range(n):
            M_new[:, i] = np.convolve(M_new[:, i], pvec)[:M_new.shape[0]]
        for k in range(i+min(M.shape[1], corr_limit, n)+1 , M_new.shape[0]):
            M_new[k, i] = 0
    print(M_new)
    print(np.matmul(M_new, M))
    #np.matmul(M_new, M)
    return M_new


def make_correction_matrix(trac_atom, formuladict, na_dict, indist_elems, autodetect, corr_limit):
    """create matrix M such that Mx=y where y is the observed isotopic distribution
    and x is the expected distribution of input labels

    formuladict: dict of element:number of atoms in molecule (e.g. {'C':2,'O':1,'H':6})
    trac_atom: element with input labeling
    indist_elems: elements with identical mass shift
    na_dict: dict of - element:expected isotopic distribution
    :TODO This function relates to issue NCT-247. Need to change the function
    in more appropriate way.
    """
    lookup_dict = {'O':['O16','O17','O18'], 'S':['S32','S33','S34'], 'Si':['Si28','Si29','Si30']}
    M = make_expected_na_matrix(formuladict.get(trac_atom, 0), na_dict[trac_atom])
    indist_elems_copy = indist_elems
    na_dict_copy = copy(na_dict)
    for e in indist_elems:
        e2 = get_isotope_element(e)
        try:
            if(lookup_dict[e2][1] in indist_elems_copy) and (lookup_dict[e2][2] in indist_elems_copy):
                # indist_elems_copy.remove(str(lookup_dict[e2][1]))
                # indist_elems_copy.remove(str(lookup_dict[e2][2]))
                # indist_elems.append(e2)
                # pos = lookup_dict[e2].index(str(e))
                # list_values = [0]*3
                # list_values[0]= na_dict_copy[e2][0]
                # list_values[pos]= na_dict_copy[e2][pos]
                na_dict_copy[lookup_dict[e2][1]] = [na_dict_copy[e2][0],na_dict_copy[e2][1],0]
                na_dict_copy[lookup_dict[e2][2]] = [0,0,na_dict_copy[e2][2]]

                #na_dict_copy[str(e)]=list_values  

            elif ((lookup_dict[e2][1] in indist_elems_copy) and (lookup_dict[e2][2] not in indist_elems_copy)) or \
                            ((lookup_dict[e2][1] not in indist_elems_copy) and (lookup_dict[e2][2] in indist_elems_copy)):
                pos = lookup_dict[e2].index(str(e))
                list_values = [0]*3
                list_values[0]= na_dict_copy[e2][0]
                list_values[pos]= na_dict_copy[e2][pos]
                na_dict_copy[str(e)]=list_values  
        except KeyError :
        	pass
        print(e2)
        print(na_dict_copy)
        M_indist = []
        if e2 in formuladict:
            e1 = ISOTOPE_NA_MASS[KEY_ELE][e]
            try:
                if autodetect:
                    print(corr_limit[e])
                    M_indist.append(add_indistinguishable_element_for_autodetect(M, formuladict[e1], na_dict_copy[e], int(corr_limit[e])))

                else:
                    M = add_indistinguishable_element(M, formuladict[e1], na_dict_copy[e])
            except:
                if autodetect:
                    print(corr_limit[e])
                    M_indist.append(add_indistinguishable_element_for_autodetect(M, formuladict[e1], na_dict_copy[e], int(corr_limit[e])))
                else:
                    M = add_indistinguishable_element(M, formuladict[e1], na_dict_copy[e2])
        if M_indist:
            i = np.eye(M.shape[1])
            for m in M_indist:
                if not np.all(m==0):
                    m = np.matmul(m, i)
                    i=m
            M = np.matmul(i, M)
    return pinv(M)


def make_all_corr_matrices(isotracers, formula_dict, na_dict, eleme_corr, autodetect, corr_limit):
    """
    This function forms correction matrix M, such that Mx=y where y is the 
    observed isotopic distribution and x is the expected distribution of input
    labels, for each indistinguishable element for a particular isotracer one by one. 
    Args:
    isotracers - list of isotracers presnt in the formula.
    formula_dict - dict of form- element:number of atoms in molecule (e.g. {'C':2,'O':1,'H':6})
    na_dict - dict of form- element:expected isotopic distribution
    eleme_corr - dict of form- isotracer_element:indistinguishable elements list(elements with identical mass shift)
                 Eg- { 'C13': ['H', 'O17'], 'N15': ['H']}

    Returns:
    corr_mats - dict of form- isotracer_element: correction matrix
    """
    corr_mats = {}
    for isotracer in isotracers:
        trac_atom = get_isotope_element(isotracer)
        try:
            indist_list = eleme_corr[trac_atom]
        except KeyError:
            indist_list = []
        corr_mats[isotracer] = make_correction_matrix(trac_atom, formula_dict, na_dict, indist_list, autodetect, corr_limit)
    return corr_mats


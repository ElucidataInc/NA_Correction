
from copy import copy

import numpy as np
from numpy.linalg import pinv
import pandas as pd
from scipy.special import binom


from corna.constants import ISOTOPE_NA_MASS, KEY_ELE
from corna.helpers import get_isotope_element

def make_expected_na_matrix(N, pvec):
    """for a single labeled element, create the matrix M
    such that Mx=y where x is the actual distribution of input labels
    and y is the expected distribution of intensities with natural abundance

    N: number of atoms of this element
    pvec: expected isotopic distribution (e.g. [0.99,0.01])"""
    ##this will work for +1 Isotopes like C13, N15, O17 but not O18
    max_label=1+(N*(len(pvec)-1))
    correction_matrix = np.zeros((max_label,N+1))
    for i in range(N+1):
        column = np.zeros(i+1)
        column[-1]=1.0
        for nb in range(N-i):
            column = np.convolve(column, pvec)
        column.resize(max_label)
        correction_matrix[:, i] = column
    return correction_matrix[:N+1]
    
def make_expected_na_matrix_3(N, pvec):
    """for a single labeled element, create the matrix M
    such that Mx=y where x is the actual distribution of input labels
    and y is the expected distribution of intensities with natural abundance

    N: number of atoms of this element
    pvec: expected isotopic distribution (e.g. [0.99,0.01])"""
    ##for +2 isotopologues
    max_label=1+(N*(len(pvec)-1))
    M_new = np.zeros((max_label,max_label))
    
    for i in range(N+1):
        for j in range(N+1):
            k=i+2*j
            if(i+j>N or k>max_label):
                break
            else:
                for m in range(max_label-k):
                    M_new[k+m][m] = M_new[k+m][m] + na_term(pvec,N,i,j)
    M_new[max_label-1][max_label-1] = 1
    return M_new


def na_term(pvec,n,i,j=0):
    if len(pvec)==2:
        coeff = binom(n,i)*(pvec[0])**(n-i)*(pvec[1])**i
    elif len(pvec)==3:
        coeff = binom(n, n-i-j)*binom(i+j,i)*pvec[0]**(n-i-j)*pvec[1]**i*pvec[2]**j
    return coeff


def add_indistinguishable_element(M,n,pvec):
    """to a matrix M formed by make_expected_na_matrix, add additional expected
    intensity corresponding to natural abundance from elements with same isotopic mass shift

    M: previous matrix formed by make_expected_na_matrix
    n: number of atoms of new element
    pvec: expected isotopic distribution of new element (e.g. [0.99,0.01])"""
    #max_label=(n*(len(pvec)-1))
    #M_new=np.zeros((M.shape[0]+max_label,M.shape[1]))
    M_new = np.zeros((M.shape[0], M.shape[1]))
    M_new[:M.shape[0],:]=M    
    for i in range(M.shape[1]):
        for j in range(n):
            M_new[:, i] = np.convolve(M_new[:, i], pvec)[:M_new.shape[0]]
    return M_new

def add_indistinguishable_element_for_autodetect_2(M, n, pvec, corr_limit):
    M_new = np.zeros((M.shape[0], M.shape[1]))
    M_new[:M.shape[0], :] = np.eye(M.shape[1])
    for i in range(M.shape[1]):
        for j in range(n):
            M_new[:, i] = np.convolve(M_new[:, i], pvec)[:M_new.shape[0]]
        for k in range(i+min(M.shape[1], corr_limit, n)+1 , M_new.shape[0]):
            M_new[k, i] = 0
    return M_new

def add_indistinguishable_element_for_autodetect_3(M, n, pvec, corr_limit_1, corr_limit_2):
    M_new = np.zeros((M.shape[0], M.shape[1]))
    for i in range(min(n, corr_limit_1)+1):
        for j in range(min(n, corr_limit_2)+1):
            k=i+2*j
            if(i+j>n or k>M.shape[1]):
                break
            else:
                for m in range(M.shape[1]-k):
                    M_new[k+m][m] = M_new[k+m][m] + na_term(pvec,n,i,j)
    return M_new


def make_correction_matrix(isotracer, formuladict, na_dict, indist_elems, autodetect, corr_limit):
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
    trac_atom = get_isotope_element(isotracer)
    if (isotracer == 'O18' or isotracer == 'S34' or isotracer == 'Si30'):
    	M = make_expected_na_matrix_3(formuladict.get(trac_atom, 0), na_dict[trac_atom])
    else:
    	M = make_expected_na_matrix(formuladict.get(trac_atom, 0), na_dict[trac_atom])
    M_indist = []

    indist_elems_copy = copy(indist_elems)
    na_dict_copy = copy(na_dict)
    for e in indist_elems:
        if e in indist_elems_copy:
            e2 = get_isotope_element(e)
            if e2 in formuladict:
                try:
                    if(lookup_dict[e2][1] in indist_elems_copy) and (lookup_dict[e2][2] in indist_elems_copy):
                        corr_limit_1 = int(corr_limit[lookup_dict[e2][1]])
                        corr_limit_2 = int(corr_limit[lookup_dict[e2][2]])
                        M_indist.append(add_indistinguishable_element_for_autodetect_3(M, formuladict[e2], na_dict_copy[e2], corr_limit_1, corr_limit_2))
                        indist_elems_copy.remove(lookup_dict[e2][1])
                        indist_elems_copy.remove(lookup_dict[e2][2])

                    elif ((lookup_dict[e2][1] in indist_elems_copy) and (lookup_dict[e2][2] not in indist_elems_copy)) or \
                                ((lookup_dict[e2][1] not in indist_elems_copy) and (lookup_dict[e2][2] in indist_elems_copy)):
                        pos = lookup_dict[e2].index(str(e))
                        list_values = [0]*3
                        list_values[0]= na_dict_copy[e2][0]
                        list_values[pos]= na_dict_copy[e2][pos]
                        na_dict_copy[str(e)]=list_values  
                        corr_limit_1 = int(corr_limit[e])
                        M_indist.append(add_indistinguishable_element_for_autodetect_2(M, formuladict[e2], na_dict_copy[e], corr_limit_1))

                    else:
                        corr_limit_1 = int(corr_limit[e])
                        M_indist.append(add_indistinguishable_element_for_autodetect_2(M, formuladict[e2], na_dict_copy[e], corr_limit_1))

                except Exception, exception_str:
                    corr_limit_1 = int(corr_limit[e])
                    M_indist.append(add_indistinguishable_element_for_autodetect_2(M, formuladict[e2], na_dict_copy[e], corr_limit_1))

    if M_indist:
        i = np.eye(M.shape[1])
        for m in M_indist:
            m = np.matmul(m, i)
            i=m
        M = np.matmul(i, M)
    return pinv(M)


def make_all_corr_matrices(isotracers, formula_dict, na_dict, autodetect, corr_limit):
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
            indist_list = corr_limit[str(trac_atom)].keys()
        except KeyError:
            indist_list = []
        corr_mats[isotracer] = make_correction_matrix(isotracer, formula_dict, na_dict, indist_list, autodetect, corr_limit[trac_atom])
    return corr_mats


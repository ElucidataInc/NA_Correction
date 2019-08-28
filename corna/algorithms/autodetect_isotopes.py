import numpy as np
from copy import copy

import constants as cs
import helpers as hl



def get_correction_limit(res, res_mw, m, instrument, mass_diff):
    if instrument=='orbitrap':
        corr_limit = np.floor(np.divide(np.multiply(1.66, np.power(m,1.5)), res*np.power(res_mw, 0.5)*mass_diff))
    elif instrument=='ft-icr':
        corr_limit = np.floor(np.divide(np.multiply(1.66, np.power(m,2)), np.multiply(res, res_mw, mass_diff)))
    return corr_limit

def get_intrument_ppm_from_res(res, res_mw, m, instrument):
    if instrument=='orbitrap':
        inst_ppm = np.multiply(10**6,(np.divide(np.multiply(1.66, np.power(m,0.5)), np.multiply(res, np.power(res_mw, 0.5)))))
    elif instrument=='ft-icr':
        inst_ppm = np.multiply(10**6,(np.divide(np.multiply(1.66, m), np.multiply(res, res_mw)))) 
    return inst_ppm

def get_ppm_required(metabolite_mass, delta_m):
    """This function calculates the ppm required to
    distinguish between the two elements in a
    metabolite.

    Args:
        formula: formula of the metabolite
        delta_m: mass diff. between the two elements

    Returns:
        required_ppm: ppm required to distinguish two elements.
    """    
    required_ppm = 1000000 * (delta_m / metabolite_mass)
    return required_ppm


def borderline_ppm_warning(ppm_user_input, required_ppm, formula, ele):
    """This function raises warning when the required ppm and
    user input ppm are nearly equal.

    Args:
        ppm_user_input: ppm of the machine used by the user
        required_ppm: ppm required to distinguish between the two elements
        formula: formula for which the elements are distinguished
        ele: element for which ppm requirements are measured

    Returns:
        boolean value: if borderline conditions are met

    """
    if (required_ppm - cs.BORDERLINE_LIMIT) <= ppm_user_input <= (required_ppm + cs.BORDERLINE_LIMIT):
        print cs.PPM_REQUIREMENT_VALIDATION + str(formula) + ':' + ele
        return True


def get_mass_diff(isotracer,element):
    """This function fetches mass difference between
    elements.

    Args:
        isotracer: labelled element
        element: indistinguishable element

    Returns:
        mass_diff: mass difference between elements

    Excepts: KeyError
        returns: None
    """
    try:
        mass_diff = cs.MASS_DIFF_DICT[isotracer][element]
        return mass_diff
    except KeyError:
        return None


def ppm_validation(ppm_user_input, required_ppm, formula, ele):
    """This function validates the ppm requirement for a particular
    element to be indistinguishable.

    Args:
        ppm_user_input: ppm of machine
        required_ppm: ppm required
        formula: formula of the metabolite
        ele: element for which validation is carried out

    Returns:
        boolean value if the requirement is met
    """
    if (borderline_ppm_warning(ppm_user_input, required_ppm, formula, ele)) \
       or (ppm_user_input > required_ppm):
       return True

       
def get_indistinguishable_ele(isotracer, res_type, formula, element, res, res_mw, instrument):
    """This function returns element which is indistinguishable for
    a particular isotracer

    Args:
        isotracer: labeled element
        element: element in the formula
        formula: formula of the metabolite
        ppm_user_input: ppm of the machine

    Returns:
        element which is indistinguishable
    """
    metabolite_mass = hl.get_mol_weight(formula)
    if res_type == 'low res':
        #not element but isotope, need to change the name
        element_symbol = hl.get_isotope_element(element)
        corr_limit = formula[element_symbol] + 100
    elif res_type == 'ultra high res':
        corr_limit = 0
    elif res_type == 'autodetect':
        mass_diff = get_mass_diff(isotracer,element)
        if mass_diff:
            required_ppm = get_ppm_required(metabolite_mass, mass_diff)
            ppm_user_input = get_intrument_ppm_from_res(res, res_mw, metabolite_mass, instrument)
            borderline_ppm_warning(ppm_user_input, required_ppm, formula, element)
            corr_limit = get_correction_limit(res, res_mw, metabolite_mass, instrument, mass_diff)
    return element, corr_limit


def get_isotope_element_list(isotracer):
    """
    This function returns list of elements present in labelled
    isotopes.
    Args:
        isotracer: list of labelled isotopes

    Returns:
        isotracer_ele_list: list of elements present
                            in isotracer.

    """
    isotracer_ele_list = []
    for label in isotracer:
        isotracer_ele = hl.get_isotope_element(label)
        isotracer_ele_list.append(isotracer_ele)
    return isotracer_ele_list


def add_isotopes_list(indis_ele_list):
    """
    This function adds isotopes of the element
    present in the formula. These isotopes are
    then tested for indistinguishable isotopes.
    Args:
        indis_ele_list: List of indistinguishable
        elements.

    Returns:
        list with isotopes of the elements.

    """
    temp_indis_ele_list = copy(indis_ele_list)
    indis_ele_list_isotopes = []
    for ele in indis_ele_list:
        if ele in cs.ISOTOPE_DICT.keys():
            temp_indis_ele_list = list((set(temp_indis_ele_list)-set([ele]))) + list(set(cs.ISOTOPE_DICT[ele]))
    return temp_indis_ele_list


def get_element_correction_dict(formula, res_type, isotracer, res, res_mw, instrument):
    """This function returns a dictionary with all isotracer elements
    as key and indistinguishable isotopes as values.

    Args:
        ppm_user_input: ppm of the machine used.
        formula: formula of the metabolites

        isotracer: labelled element which is to be corrected

    Returns:
        element_correction_dict: element correction dictionary.
    """

    element_correction_dict = {}
    formula_dict = hl.parse_formula(formula)
    ele_list = formula_dict.keys()
    isotracer_list = get_isotope_element_list(isotracer)
    isotope_ele = get_isotope_element_list(cs.MASS_DIFF_DICT.keys())

    ele_list_without_isotracer = set(ele_list) - set(isotracer_list)

    unlabeled_ele = {}
    
    if res_type == 'low res':
        if not len(isotracer) == 1:
            raise ValueError('Multiple Isotracers not accepted for Low Res Instruments')

    elif res_type == 'ultra high res':
        if not len(isotracer) == 1:
            #removing all dictionary except for the first key
            #to avoid repeated multiplication
            for i in range(1,len(isotracer)):
                if isotracer[i][0] in ele_list:
                    element_correction_dict[isotracer[i][0]] = {}
            isotracer = [isotracer[0]]


    for isotope in isotracer:
        correction_limit_dict = {}

        #this isotope[0] needs to change, hard to read code
        #also errornous if isotracer element has two letter like Si
        if isotope[0] in ele_list:
            ##restricting correction elements to common isotopes CHNOPS, 
            ##should be given as a warning somewhere or specified in documentation
            indis_ele_list = list(ele_list_without_isotracer.intersection(set(isotope_ele)))
            indis_ele_list = add_isotopes_list(indis_ele_list)
            get_ele = lambda iso: get_indistinguishable_ele(isotope, res_type, formula_dict, iso, res, res_mw, instrument)
            indis_element = map(get_ele, indis_ele_list)
            indis_element = filter(None, indis_element)
            indis_element_list = []
            for ele, num in indis_element:
                if res_type == 'autodetect':
                    try:
                        if (unlabeled_ele[ele][1] == 0 and num == 0):
                            pass
                        elif (unlabeled_ele[ele][1] == 0 and num != 0):
                            old_iso = unlabeled_ele[ele][0]
                            element_correction_dict[old_iso].pop(ele, None)
                            correction_limit_dict[ele] = num
                            unlabeled_ele[ele] = [isotope[0], num]
                        elif (unlabeled_ele[ele][1] != 0 and num==0):
                            pass
                        elif (unlabeled_ele[ele][1] != 0 and num != 0):
                            old_iso = unlabeled_ele[ele][0]
                            element_correction_dict[old_iso].pop(ele, None)
                            #removing common indistinguishable elements
                            #user can correct on their own if they want
                            #package doesn't handle that case
                            #raise warning letting user know about this scenario
                    except:
                        unlabeled_ele[ele] = [isotope[0], num]
                        correction_limit_dict[ele] = num
                else:
                    correction_limit_dict[ele] = num
            element_correction_dict[isotope[0]] = correction_limit_dict
    return element_correction_dict


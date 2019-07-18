import numpy as np

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
        print cs.PPM_REQUIREMENT_VALIDATION + formula + ':' + ele
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

       
def get_indistinguishable_ele(isotracer, formula, element, res, res_mw, instrument):
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
    temp_indis_ele_list = indis_ele_list
    indis_ele_list_isotopes = []
    for ele in indis_ele_list:
        if ele in cs.ISOTOPE_DICT.keys():
            indis_ele_list_isotopes = list((set(temp_indis_ele_list)-set([ele]))) + list(set(cs.ISOTOPE_DICT[ele]))
            temp_indis_ele_list = indis_ele_list_isotopes
    return indis_ele_list_isotopes


def get_element_correction_dict(formula, isotracer, res, res_mw, instrument):
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
    correction_limit_dict = {}
    ele_list = (hl.parse_formula(formula)).keys()
    isotracer_list = get_isotope_element_list(isotracer)
    isotope_ele = get_isotope_element_list(cs.MASS_DIFF_DICT.keys())
    ele_list_without_isotracer = set(ele_list) - set(isotracer_list)
    for isotope in isotracer:
        if isotope[0] in ele_list:
            indis_ele_list = list(ele_list_without_isotracer.intersection(set(isotope_ele)))
            indis_ele_list = add_isotopes_list(indis_ele_list)
            get_ele = lambda iso: get_indistinguishable_ele(isotope, formula, iso, res, res_mw, instrument)
            indis_element = map(get_ele, indis_ele_list)
            indis_element = filter(None, indis_element)
            indis_element_list = []
            for ele, num in indis_element:
                correction_limit_dict[ele] = num
                indis_element_list.append(ele)
            element_correction_dict[isotope[0]] = indis_element_list
    return element_correction_dict, correction_limit_dict


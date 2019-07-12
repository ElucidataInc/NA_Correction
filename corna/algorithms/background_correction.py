import pandas as pd
import math
import numpy as np
from scipy.misc import comb
from functools import partial

from corna import helpers
from corna import data_model
from corna.helpers import get_isotope_na
from corna import constants
from corna.isotopomer import bulk_insert_data_to_fragment, Infopacket
from corna.inputs.column_conventions import multiquant
from corna.inputs import multiquant_parser

def background_noise(unlabel_intensity, na, parent_atoms, parent_label,
                                             daughter_atoms, daughter_label):
    """
    This function returns noise for a particular cohort group calculated through the
    intensity of the unlabeled fragment(C12 PARENT) of that cohort group.

    Parameters
    ----------
    unlabel_intensity : int
        intensity of M+0 isotopologue
    na : int
        na value of the isotope
    parent_atoms : int
        total no. of atoms of that element in parent fragment.
    parent_label : int
        No. of labeled atoms of that element in parent fragment.
    daughter_atoms : int
        Total no. of atoms of that element in daughter fragment.
    daughter_label : int
        No. of labeled atoms of that element in parent fragment.

    Returns
    -------
    noise : int    

    """
    noise = unlabel_intensity * math.pow(na, parent_label)\
        * comb(parent_atoms - daughter_atoms, parent_label - daughter_label)\
        * comb(daughter_atoms, daughter_label)
    return noise


def background_subtraction(input_intensity, noise):
    """
    This function returns background noise of each sample.
    If background noise is negative it returns zero.

    Parameters
    ----------
    input_intensity : int
        intensity of sample
    noise : int 
        noise of cohort group    
    """
    intensity = input_intensity - noise
    return intensity
    # if intensity > 0:
    #     return intensity
    # else:
    #     return 0


def background(list_of_replicates, input_fragment_value, unlabeled_fragment_value, isotope_dict):
    parent_frag, daughter_frag = input_fragment_value.frag
    iso_elem = helpers.get_isotope_element(parent_frag.isotracer)
    parent_label = parent_frag.get_num_labeled_atoms_isotope(
        parent_frag.isotracer)
    parent_atoms = parent_frag.number_of_atoms(iso_elem)
    na = helpers.get_isotope_na(parent_frag.isotracer, isotope_dict)
    daughter_atoms = daughter_frag.number_of_atoms(iso_elem)
    daughter_label = daughter_frag.get_num_labeled_atoms_isotope(
        parent_frag.isotracer)
    replicate_value = {}
    for replicate_group in list_of_replicates:

        background_list = []
        for each_replicate in replicate_group:
            noise = background_noise(unlabeled_fragment_value.data[each_replicate], na, parent_atoms,
                                     parent_label, daughter_atoms, daughter_label)
            try:
                background = background_subtraction(input_fragment_value.data[each_replicate], noise)
            except:
                background = background_subtraction(0, noise)
            background_list.append(background)
        background_value = max(background_list)
        for each_replicate in replicate_group:
            replicate_value[each_replicate] = background_value
    return replicate_value


def bacground_correction(replicates, sample_background, sample_data, decimals):
    corrected_sample_data = {}
    for key, value in sample_data.iteritems():
        background_value = replicates[sample_background[key]]
        new_value = np.around(value - background_value, decimals)
        corrected_sample_data[key] = new_value
    return corrected_sample_data


def bulk_background_correction(fragment_dict, list_of_replicates, sample_background, isotope_dict, decimals):
    input_fragments = []
    unlabeled_fragment = []
    corrected_fragments_dict = {}
    for key, value in fragment_dict.iteritems():
        if value.unlabeled:
            unlabeled_fragment.append((key, value))
            input_fragments.append((key, value))
        else:
            input_fragments.append((key, value))
    try:
        assert len(unlabeled_fragment) == 1
    except AssertionError:
        raise AssertionError('The input should contain atleast and only one unlabeled fragment data'
                             'Please check metadata or raw data files')
    for input_fragment in input_fragments:
        replicate_value = background(list_of_replicates, input_fragment[
            1], unlabeled_fragment[0][1], isotope_dict)
        #corrected_sample_data = bacground_correction(
        #    replicate_value, sample_background, input_fragment[1].data, decimals)
        sample_data = input_fragment[1].data
        corrected_sample_data = {}
        for key, value in sample_data.iteritems():
            background_value = replicate_value[sample_background[key]]
            new_value = np.around(value - background_value, decimals)
            corrected_sample_data[key] = new_value

            corrected_fragments_dict[input_fragment[0]] = Infopacket(input_fragment[1].frag, corrected_sample_data,
                                                       input_fragment[1].unlabeled, input_fragment[1].name)
    return corrected_fragments_dict


# def background_correction(metabolite_frag_dict, list_of_replicates, sample_background, isotope_dict=constants.ISOTOPE_NA_MASS, decimals=0):
#     preprocessed_output_dict = {}
#     for metabolite, fragments_dict in metabolite_frag_dict.iteritems():
#         preprocessed_output_dict[metabolite] = bulk_background_correction(fragments_dict,
#                                                                           list_of_replicates,
#                                                                           sample_background, isotope_dict, decimals)

#     return preprocessed_output_dict

def background_correction(msms_df, list_of_replicates, sample_background, isotope_dict=constants.ISOTOPE_NA_MASS, decimals=0):
    metabolite_frag_dict = multiquant_parser.mq_df_to_fragmentdict(msms_df)
    preprocessed_output_dict = {}
    for metabolite, fragments_dict in metabolite_frag_dict.iteritems():
        preprocessed_output_dict[metabolite] = bulk_background_correction(fragments_dict,
                                                                          list_of_replicates,
                                                                          sample_background, isotope_dict, decimals)

    return preprocessed_output_dict


import pandas as pd
import math
import numpy as np
from scipy.misc import comb

from corna.helpers import get_isotope_na
from corna import constants as const
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
    if intensity > 0:
        return intensity
    else:
        return 0


def background_correction(msms_df, list_of_replicates, isotope_dict=const.ISOTOPE_NA_MASS):
    """
    This function corrects intensity value errors present due to background noise.
    
    Parameters
    ----------
    msms_df:
         Dataframe which contains intensities to be corrected.
    list_of_replicates:
         List consisting all the sample group for each cohort.

    Returns
    -------
        output_df: Background corrected dataframe
    """
    final_df, isotracer = multiquant_parser.add_mass_and_no_of_atoms_info_frm_label(msms_df)
    output_df = pd.DataFrame()

    na = get_isotope_na(isotracer[0], isotope_dict)

    final_df[const.UNLABELED] = 'False'
    final_df.loc[final_df[const.PARENT_NUM_LABELED_ATOMS] == 0, const.UNLABELED] = 'True'

    for metab in final_df.Name.unique():
        metab_df = final_df[final_df[const.NAME_COL] == metab]
        unlabel_isotope_df = metab_df[metab_df[const.UNLABELED] == 'True']

        for frag in metab_df[multiquant.MQ_FRAGMENT].unique():
            frag_df = metab_df[metab_df[multiquant.MQ_FRAGMENT] == frag]
            for replicate_group in list_of_replicates:
                background_list = []
                for each_replicate in replicate_group:
                    unlabel_intensity = unlabel_isotope_df.loc[unlabel_isotope_df[const.SAMPLE_COL] == each_replicate,
                                                                         const.INTENSITY_COL].iloc[0]
                    noise = background_noise(unlabel_intensity, na, frag_df[const.PARENT_NUM_ATOMS].unique()[0],
                                                    frag_df[const.PARENT_NUM_LABELED_ATOMS].unique()[0],
                                                    frag_df[const.DAUGHTER_NUM_ATOMS].unique()[0],
                                                    frag_df[const.DAUGHTER_NUM_LABELED_ATOMS].unique()[0])
                    background = background_subtraction(frag_df.loc[frag_df[const.SAMPLE_COL] == each_replicate, 
                                                                    const.INTENSITY_COL].iloc[0], noise)
                    background_list.append(background)
                background_value = max(background_list)
                for each_replicate in replicate_group:
                    frag_df.loc[frag_df[multiquant.BACKGROUND] == each_replicate, 'replicate_value'] = background_value
            output_df= output_df.append(frag_df)

    output_df[const.BACKGROUND_CORRECTED] = output_df[const.INTENSITY_COL] - output_df['replicate_value']
    output_df.drop('replicate_value', axis=1, inplace=True)
    return output_df


def replace_negatives_background(df):
    """
    This function replaces negative intensity values to zero.
    """
    df[const.BACKGROUND_WITH_ZERO] = df[const. BACKGROUND_CORRECTED].clip(lower=0)
    return df

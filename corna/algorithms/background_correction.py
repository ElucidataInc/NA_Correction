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
    return intensity
    #if intensity > 0:
    #    return intensity
    #else:
    #    return 0


def background_correction(msms_df, list_of_replicates, isotope_dict=const.ISOTOPE_NA_MASS):
    """
    This function corrects intensity value errors present due to background noise.
    
    Parameters
    ----------
    msms_df:
         Dataframe which contains intensities to be corrected.
    list_of_replicates:
         List consisting all the sample group for each cohort.
    isotope_dict:
         Dictionary containing the natural abundance values of each element and its isotope.
         Eg: {'C12': 0.9889, 'C13':0.0111, 'N14':0.9964,'N15':0.0036}
    Returns
    -------
        output_df: Background corrected dataframe
    """
    final_df, isotracer = multiquant_parser.add_mass_and_no_of_atoms_info_frm_label(msms_df)
    #print(final_df.columns)
    #print(isotracer)
    output_df = pd.DataFrame()

    na = get_isotope_na(isotracer[0], isotope_dict)
    #print(na)

    final_df[const.UNLABELED] = 'False'
    final_df.loc[final_df[const.PARENT_NUM_LABELED_ATOMS] == 0, const.UNLABELED] = 'True'
    #rint(final_df.head())
    final_df.to_csv('final_df.csv')

    for metab in final_df.Name.unique():
        #print(metab)
        metab_df = final_df[final_df[const.NAME_COL] == metab]
        unlabel_isotope_df = metab_df[metab_df[const.UNLABELED] == 'True']
        metab_df.to_csv('metab_df.csv')
        unlabel_isotope_df.to_csv('unlabel_isotope_df.csv')

        for frag in metab_df[multiquant.MQ_FRAGMENT].unique():
            #print(frag)
            frag_df = metab_df[metab_df[multiquant.MQ_FRAGMENT] == frag]
            #print(frag_df.head(2))
            frag_df.to_csv('frag_df.csv')
            for replicate_group in list_of_replicates:
                background_list = []
                for each_replicate in replicate_group:
                    try:
                        unlabel_intensity = unlabel_isotope_df.loc[unlabel_isotope_df[const.SAMPLE_COL] == each_replicate,
                                                                         const.INTENSITY_COL].iloc[0]                                           
                    except:
                        unlabel_intensity = 0
                    #print('unlabel intensity') 
                    #print(each_replicate, unlabel_intensity) 
                    noise = background_noise(unlabel_intensity, na, frag_df[const.PARENT_NUM_ATOMS].unique()[0],
                                                    frag_df[const.PARENT_NUM_LABELED_ATOMS].unique()[0],
                                                    frag_df[const.DAUGHTER_NUM_ATOMS].unique()[0],
                                                    frag_df[const.DAUGHTER_NUM_LABELED_ATOMS].unique()[0])
                    try:
                        intensity = frag_df.loc[frag_df[const.SAMPLE_COL] == each_replicate, 
                                                                    const.INTENSITY_COL].iloc[0]
                    except:
                        intensity = 0
                    background = background_subtraction(intensity , noise)
                    background_list.append(background)
                background_value = max(background_list)
                for each_replicate in replicate_group:
                    frag_df.loc[frag_df[multiquant.BACKGROUND] == each_replicate, 'replicate_value'] = background_value
            output_df= output_df.append(frag_df)

    output_df[const.BACKGROUND_CORRECTED] = output_df[const.INTENSITY_COL] - output_df['replicate_value']
    output_df.drop('replicate_value', axis=1, inplace=True)
    return output_df



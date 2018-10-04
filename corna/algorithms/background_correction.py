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
    isotope_dict:
         Dictionary containing the natural abundance values of each element and its isotope.
         Eg: {'C12': 0.9889, 'C13':0.0111, 'N14':0.9964,'N15':0.0036}
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
                    try:
                        unlabel_intensity = unlabel_isotope_df.loc[unlabel_isotope_df[const.SAMPLE_COL] == each_replicate,
                                                                         const.INTENSITY_COL].iloc[0]
                    except:
                        unlabel_intensity = 0
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


msms_df = pd.read_csv('/home/priyanka/Downloads/raaaaaaaaaaaaisa.csv')
print msms_df
list_of_replicates = [[ 'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (10) | Sample Number=61 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (11) | Sample Number=67 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (12) | Sample Number=73 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (13) | Sample Number=79 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (8) | Sample Number=49 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (9) | Sample Number=55 '], [ 'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (2) | Sample Number=13 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (3) | Sample Number=19 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (4) | Sample Number=25 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (5) | Sample Number=31 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (6) | Sample Number=37 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc B6 0 1-20dil (7) | Sample Number=43 '], [ 'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (10) | Sample Number=133 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (11) | Sample Number=139 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (12) | Sample Number=145 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (13) | Sample Number=151 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (8) | Sample Number=121 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (9) | Sample Number=127 '], [ 'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (2) | Sample Number=85 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (3) | Sample Number=91 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (4) | Sample Number=97 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (5) | Sample Number=103 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (6) | Sample Number=109 ',
       'Wisconsin_rbcplates_10aug18_11sep18-rbc PWK 0 1-20dil (7) | Sample Number=115 ']]
final_df = background_correction(msms_df, list_of_replicates, isotope_dict=const.ISOTOPE_NA_MASS)

final_df.to_csv('/home/priyanka/Downloads/raisa2.csv')
import pandas as pd
import math
import numpy as np
from scipy.misc import comb

from corna.helpers import get_isotope_na
from corna import constants as const
from corna.inputs.column_conventions import multiquant
from corna.inputs import multiquant_parser 

def background_noise(unlabel_intensity, na, parent_atoms, parent_label, daughter_atoms, daughter_label):
    """
    This function returns noise due to each sample for a particular cohort.
    """
    noise = unlabel_intensity * math.pow(na, parent_label)\
        * comb(parent_atoms - daughter_atoms, parent_label - daughter_label)\
        * comb(daughter_atoms, daughter_label)
    return noise

def backround_subtraction(input_intensity, noise):
    intensity = input_intensity - noise
    return intensity

def background_correction(msms_df, list_of_replicates):
    """
    This function corrects intensity value errors present due to background noise.
    Args:
        msms_df: Dataframe which contains intensities to be corrected.
        list_of_replicates: List consisting all the sample group for each cohort.

    Returns:
        output_df: Background corrected dataframe
    """
    final_df, isotracer= multiquant_parser.add_info_to_df(msms_df)
    output_df= pd.DataFrame()

    na= get_isotope_na(isotracer[0])

    final_df[const.UNLABELED]= 'False'
    final_df.loc[final_df[const.PARENT_NUM_LABELED_ATOMS] == 0, const.UNLABELED] = 'True'
    
    for metab in final_df.Name.unique():
        metab_df= final_df[final_df[const.NAME_COL]== metab]
        calcu_df= metab_df[metab_df[const.UNLABELED]== 'True']

        for frag in metab_df[multiquant.MQ_FRAGMENT].unique():
            frag_df= metab_df[metab_df[multiquant.MQ_FRAGMENT]== frag]
            for replicate_group in list_of_replicates:
                background_list = []
                for each_replicate in replicate_group:
                    unlabel_intensity= calcu_df.loc[calcu_df[const.SAMPLE_COL] == each_replicate, const.INTENSITY_COL].iloc[0]
                    noise = background_noise(unlabel_intensity, na, frag_df[const.PARENT_NUM_ATOMS].unique()[0], frag_df[const.PARENT_NUM_LABELED_ATOMS].unique()[0],
                                     frag_df[const.DAUGHTER_NUM_ATOMS].unique()[0], frag_df[const.DAUGHTER_NUM_LABELED_ATOMS].unique()[0])
                    background = backround_subtraction(frag_df.loc[frag_df[const.SAMPLE_COL] == each_replicate, const.INTENSITY_COL].iloc[0], noise)
                    background_list.append(background)
                background_value = max(background_list)
                for each_replicate in replicate_group:
                    frag_df.loc[frag_df[multiquant.BACKGROUND] == each_replicate, 'replicate_value'] = background_value
            output_df= output_df.append(frag_df)

    output_df[const.BACKGROUND_CORRECTED]= output_df[const.INTENSITY_COL]- output_df['replicate_value']
    output_df= replace_negatives_background(output_df)
    return output_df

def replace_negatives_background(df):
    """
    This function replaces negative intensity values to zero.
    """
    df[const.BACKGROUND_WITH_ZERO]= df[const. BACKGROUND_CORRECTED].clip(lower=0)
    return df  

        





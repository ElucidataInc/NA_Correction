from collections import namedtuple
import pandas as pd
import numpy as np
import re

import corna.constants as const
from corna.helpers import get_isotope_na
from corna.inputs.column_conventions import multiquant 
from corna.inputs import multiquant_parser


def na_correction_mimosa(msms_df, isBackground):
    """
    This function performs na correction on the input data frame for LCMS/MS file.
    Args:
        msms_df: Dataframe which contains intensities to be corrected.

    Returns:
        output_df: na corrected dataframe
        
    """
    REQUIRED_COL= [multiquant.FORMULA, multiquant.LABEL, multiquant.NAME, multiquant.SAMPLE, multiquant.COHORT,
                     multiquant.MQ_FRAGMENT, multiquant.INTENSITY, multiquant.PARENT_FORM,const.NA_CORRECTED_COL]
                     
    if isBackground:
        final_df= msms_df
        isotracer= msms_df[multiquant.ISOTRACER].unique()
        intensity_col= const.BACKGROUND_WITH_ZERO
    else:        
        final_df, isotracer= multiquant_parser.add_info_to_df(msms_df)
        intensity_col= const.INTENSITY_COL

    final_df[const.NA_CORRECTED_COL]=0.0
    output_df= pd.DataFrame()
    metab_dict={}

    na= hlp.get_isotope_na(isotracer[0])

    final_df['A']=(1 + na * (final_df[const.PARENT_NUM_ATOMS]-final_df[const.PARENT_NUM_LABELED_ATOMS]))
    final_df['B']= na * ((final_df[const.PARENT_NUM_ATOMS]-final_df[const.DAUGHTER_NUM_ATOMS]) -\
                         (final_df[const.PARENT_NUM_LABELED_ATOMS]-final_df[const.DAUGHTER_NUM_LABELED_ATOMS]-1))
    final_df['C']=  na * (final_df[const.DAUGHTER_NUM_ATOMS]-final_df[const.DAUGHTER_NUM_LABELED_ATOMS]+1)

    final_df.drop([const.PARENT_MASS_MOL, const.DAUGHTER_MASS_MOL, const.PARENT_NUM_ATOMS,
                const.DAUGHTER_NUM_ATOMS, const.DAUGHTER_NUM_LABELED_ATOMS, const.PARENT_NUM_LABELED_ATOMS], axis=1, inplace=True)

    for samp in final_df.Sample.unique():
        """
        Create metabolite dictionary of the form:
            {'SAMPLE 2_10':{
                (191, 111): 2345.75, (192, 111):5644.847
                }
            }
        """
        metab_df = final_df[final_df[multiquant.SAMPLE]==samp]
        frag_dict={}
        for index, row in metab_df.iterrows():
            frag_dict[(row[const.PARENT_MASS_ISO],row[const.DAUGHTER_MASS_ISO])]=row[intensity_col]
    
        metab_dict[samp]= frag_dict

    for samp in final_df.Sample.unique():
        metab_df = final_df[final_df[multiquant.SAMPLE]==samp]
        frag= metab_dict[samp]  
     
        for index, row in metab_df.iterrows():
            m_n= row[const.DAUGHTER_MASS_ISO]
            m_1_n= row[const.PARENT_MASS_ISO]-1
            m_n_1= row[const.DAUGHTER_MASS_ISO]-1
            intensity_m_n= row[intensity_col]
            try:
                intensity_m_1_n= frag[m_1_n, m_n]
            except KeyError:
                intensity_m_1_n=0
            try:
                intensity_m_1_n_1= frag[m_1_n, m_n_1]
            except KeyError:
                intensity_m_1_n_1= 0
        
            corrected= intensity_m_n * row['A']  - intensity_m_1_n * row['B'] -\
                                                intensity_m_1_n_1 * row['C']
            metab_df.set_value(index=index, col=const.NA_CORRECTED_COL, value= corrected)

        output_df=output_df.append(metab_df) 
    output_df= output_df.filter(REQUIRED_COL)

    return output_df



  


        






    
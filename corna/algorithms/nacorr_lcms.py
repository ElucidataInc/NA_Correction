"""
This module is a wrapper around algorithms.py. It calls different functions from algorithms.py
and performs na correction (na_correction function). The output is given in the form of dataframe 
which can further be used in post processing function etc.
"""

import numpy as np
import pandas as pd

import matrix_calc as algo
import corna.inputs.maven_parser as parser
from corna.autodetect_isotopes import get_element_correction_dict
from corna.constants import INTENSITY_COL, ISOTOPE_NA_MASS, KEY_NA
from corna.helpers import get_isotope_element, first_sub_second, parse_formula, chemformula_schema

def eleme_corr_invalid_entry(iso_tracers, eleme_corr):
    """
    This function raises an error if the user inputs incorrect values in  eleme_corr dictionary.
    The indistinguishable element specified in the eleme_corr dictionary cannot be an isotracer
    element. This is logically incorrect input, hence raises an error.
    """
    for key, value in eleme_corr.iteritems():
        for isotracer in iso_tracers:
            el = get_isotope_element(isotracer)
            if el in value:
                raise KeyError('An iso tracer cannot'
                               ' be an Indistinguishable element (' + el +
                               ') , invalid input in eleme_corr dictionary')


def correct_df_for_multiplication(isotracers, lab_samp_df, corr_mats):
    """
    This function arranges the dataframe according to each isotope by grouping it on the basis of that isotope column 
    and then pass it as input to perform multiplication with correction matrix.
    """
    curr_df = lab_samp_df
    if len(isotracers) == 1:
        curr_df = multiplying_df_with_matrix(isotracers[0], corr_mats[isotracers[0]], curr_df)
    else:
        for isotracer in isotracers:
            index_cols = first_sub_second(isotracers, [isotracer])
            group_by_index_cols = curr_df.groupby(level = index_cols)
            L=[]
            keys=[]
            for group_no, group in group_by_index_cols:
                group.index = group.index.get_level_values(isotracer)
                corr_df = multiplying_df_with_matrix(isotracer, corr_mats[isotracer], group)
                L.append(corr_df)
                keys.append(group_no)
            curr_df = pd.concat(L, keys=keys, names=index_cols)
            curr_df.index = curr_df.index.reorder_levels(lab_samp_df.index.names)

    return curr_df
    

def multiplying_df_with_matrix(isotracer, corr_mat_for_isotracer, curr_df):
    """This function takes the correction matrix for given isotracer and multiplies
    it with the sample values of the dataframe to give corrected sample values
    Example:
        C13 Sample1 Sample2
        0   0.21    0.98
        1   0.34    0.11
        multiplied with correction matrix [0.99 0]
                                          [0.01 1]
        will give output

        C13 Sample1 Sample2
         0  0.2079  0.9702
	     1  0.3421  0.1198
    """
    num_rows, num_cols = corr_mat_for_isotracer.shape
    curr_df = curr_df.reindex(np.arange(num_cols)).fillna(0)
    corr_data = np.matmul(corr_mat_for_isotracer, curr_df.values)
    corr_df = pd.DataFrame(index=pd.index.np.arange(num_rows),columns=curr_df.columns, data=corr_data)
    corr_df.index.name = isotracer
    return corr_df


def na_correction(merged_df, iso_tracers, ppm_input_user, na_dict, eleme_corr, intensity_col=INTENSITY_COL,autodetect=False):


    """
    This function performs na correction on the input data frame for LCMS file. 
    Args:
        merged_df: data frame which contains intensities which are to be corrected
        iso_tracers: list of labeled elements. eg ['C13', 'N15']
        ppm_input_user: ppm resolution of the machine required when user
                        selects autodetect = True, else a blank parameter
                        can be given.
        na_dict: dictionary with natural abundance values of the elements.
        eleme_corr: if user selects autodetect=False, they can give a standard
                    dict of indistinguishable elements for correction.
                    eg - {'C13':['H','O']}
        autodetect:It takes boolean value for auto detection. By default it is False.

    Returns:
        joined: na corrected dataframe
        eleme_corr_dict : dictionary od indistinguishable isotopes used for correction
    """
    original_df= parser.save_original_df(merged_df, iso_tracers)
    
    sample_list=merged_df.Sample.unique()
    required_col=np.append(sample_list, iso_tracers)
    final_df=pd.DataFrame()
    eleme_corr_dict = {}

    #convert input datframe from long to wide format
    merged_df=merged_df.pivot_table(index=['Name','Formula','Label'], columns='Sample', values='Intensity')
    merged_df =merged_df.rename_axis(None, axis=1).reset_index()
    
    std_label_df = parser.convert_labels_to_std(merged_df, iso_tracers)   
     
    if autodetect:
        for metab in std_label_df.Name.unique():
            
            required_df, formula, formula_dict = parser.process_corrected_df_for_metab(std_label_df, metab, 
                                                                    iso_tracers,required_col)
            auto_eleme_corr = get_element_correction_dict(ppm_input_user, formula[0] ,iso_tracers)
            eleme_corr_dict[metab] = auto_eleme_corr
            corr_mats = algo.make_all_corr_matrices(iso_tracers, formula_dict, na_dict, eleme_corr)
            df_corr_C_N = correct_df_for_multiplication(iso_tracers, required_df, corr_mats)
            info_df= parser.add_info_to_final_df(df_corr_C_N, metab, formula[0], iso_tracers)
            final_df=final_df.append(info_df)
    
    else:
        eleme_corr_invalid_entry(iso_tracers, eleme_corr)
        
        for metab in std_label_df.Name.unique():
            eleme_corr_dict[metab] = eleme_corr  
            required_df, formula, formula_dict = parser.process_corrected_df_for_metab(std_label_df, metab,
                                                                     iso_tracers, required_col)
            corr_mats = algo.make_all_corr_matrices(iso_tracers, formula_dict, na_dict, eleme_corr)
            df_corr_C_N = correct_df_for_multiplication(iso_tracers, required_df, corr_mats)
            info_df= parser.add_info_to_final_df(df_corr_C_N, metab, formula[0], iso_tracers)
            final_df=final_df.append(info_df)
        
    #convert na corrected datframe back from wide format to long format        
    df_long = pd.melt(final_df, id_vars=['Name', 'Formula', 'Label'])
    df_long.rename(columns={'variable': 'Sample', 'value':'NA Corrected'},inplace=True)

    #merge original_df with na corrected df to get specific column information
    joined= pd.merge(df_long, original_df, left_on=['Name', 'Formula', 'Label', 'Sample'], right_on=['Name', 'Formula','Label', 'Sample'], how='left')        
    joined= joined.fillna(0)
    joined.loc[joined.Original_label == 0, 'Original_label'] = joined.Label
    joined.drop('Label', axis=1, inplace=True)
    joined.rename(index= str, columns={'Original_label':'Label'}, inplace=True)    
    return joined, eleme_corr_dict

    

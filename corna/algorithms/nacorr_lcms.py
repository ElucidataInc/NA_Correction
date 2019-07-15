"""
This module is a wrapper around algorithms.py. It calls different functions from algorithms.py
and performs na correction (na_correction function). The output is given in the form of dataframe 
which can further be used in post processing function etc.
"""
from copy import copy

import numpy as np
import pandas as pd

import matrix_calc as algo
import corna.inputs.maven_parser as parser
from corna.autodetect_isotopes import get_element_correction_dict
from corna import constants as cons
from corna.helpers import get_isotope_element, first_sub_second, parse_formula, chemformula_schema, get_na_value_dict

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


def get_correct_df_by_multiplication(isotracers, preprocessed_df, corr_mats):
    """
    This function arranges the dataframe according to each isotope by grouping it on the basis of that isotope column 
    and then pass it as input to perform multiplication with correction matrix.
    Finally NACorrected df is obtained.
    """
    curr_df = preprocessed_df
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
            curr_df.index = curr_df.index.reorder_levels(preprocessed_df.index.names)

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
    corr_df = pd.DataFrame(data=corr_data, index=pd.index.np.arange(num_rows),columns=curr_df.columns)
    corr_df.index.name = isotracer
    return corr_df

def perform_nacorrection_metab(df, metab, iso_tracers, required_col, na_dict, eleme_corr,final_df, autodetect, corr_limit):
    """
    This function performs na correcion for each metabolite one by one, adds required info
    back to the dataframe and then returns the NA corrected dataframe.
    To perform NA correction, first formula dictionary of th metabolite is calculated
    ,which is the number of each element in the compound.Then correction matrix is calculated 
    based on the formula and isotracer. The values of Dataframe is then multiplied to the corresponding 
    values of the dataframe to get nacorrected values.These nacorrected values are then appended in
    the dataframe.
    Args:
        df: data frame of one metabolite which contains intensities which are to be corrected
        metab: metabolite name whose dataframe is passed
        iso_tracers: list of labeled elements. eg ['C13', 'N15']
        required_col: list of column names required for na correction calculations.
        na_dict: dictionary with natural abundance values of the elements.
        eleme_corr: if user selects autodetect=False, they can give a standard
                    dict of indistinguishable elements for correction.
                    eg - {'C13':['H','O']}
        final_df: NA corrected df of other metabolites which needs to be appended to the calculated df.

    Returns:
        final_df: NA corrected dataframe of a single metabolite                
    
    """
    required_df, formula, formula_dict = parser.filter_required_col_and_get_formula_dict(df, metab,
                                                                     iso_tracers, required_col)
    corr_mats = algo.make_all_corr_matrices(iso_tracers, formula_dict, na_dict, eleme_corr, autodetect, corr_limit)
    corrected_df = get_correct_df_by_multiplication(iso_tracers, required_df, corr_mats)
    info_df= parser.add_name_formula_label_col(corrected_df, metab, formula[0], iso_tracers, eleme_corr)
    final_df=final_df.append(info_df)
    return final_df


def na_correction(merged_df, iso_tracers, eleme_corr, na_dict=get_na_value_dict(), autodetect=False, res=None, res_mw=None, instrument=None):
    """
    This function performs na correction on the input data frame for LCMS file. 
    This function is a wrapper around perform_nacorrection_metab function. It preprocesses 
    the input to get the desired input format of dataframe for each metabolite one by one and 
    then pass the processed dataframe to the perform_nacorrection_metab to get nacorrected 
    dataframe of a single metabolite. And keeps on appending this nacorrected dataframe to get 
    the final joined dataframe. 
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
        autodetect:It takes boolean value for auto detection of Indistinguishable Isotopes. 
                   By default it is False.

    Returns:
        joined: na corrected dataframe
        eleme_corr_dict : dictionary od indistinguishable isotopes used for correction
    """
    original_df= parser.save_original_label_and_processed_label(merged_df, iso_tracers)
    
    sample_list=merged_df.Sample.unique()
    required_col=np.append(sample_list, iso_tracers)
    final_df=pd.DataFrame()
    eleme_corr_dict = {}

    #convert input datframe from long to wide format
    merged_df=merged_df.pivot_table(index=[cons.NAME_COL, cons.FORMULA_COL, cons.LABEL_COL],
                                                                     columns=cons.SAMPLE_COL, values= cons.INTENSITY_COL)
    merged_df =merged_df.rename_axis(None, axis=1).reset_index()
    
    std_label_df = parser.get_isotope_columns_frm_label_col(merged_df, iso_tracers)   
    if autodetect:
        for metab in std_label_df.Name.unique():
            formula= std_label_df[std_label_df[cons.NAME_COL]== metab].Formula.unique()
            auto_eleme_corr, corr_limit = get_element_correction_dict(formula[0] ,iso_tracers, res, res_mw, instrument)
            eleme_corr_dict[metab] = auto_eleme_corr
            final_df= perform_nacorrection_metab(std_label_df, metab, iso_tracers, required_col, na_dict,
                                                     auto_eleme_corr, final_df, autodetect=True, corr_limit=corr_limit)            
    
    else:
        eleme_corr_invalid_entry(iso_tracers, eleme_corr)
        
        for metab in std_label_df.Name.unique():
            eleme_corr_dict[metab] = eleme_corr 
            final_df= perform_nacorrection_metab(std_label_df, metab, iso_tracers, required_col, na_dict,
                                                     eleme_corr, final_df, autodetect=False, corr_limit=None) 
            
    #convert na corrected datframe back from wide format to long format        
    df_long = pd.melt(final_df, id_vars=[cons.NAME_COL, cons.FORMULA_COL, cons.LABEL_COL, cons.INDIS_ISOTOPE_COL])
    df_long.rename(columns={'variable': cons.SAMPLE_COL, 'value':cons.NA_CORRECTED_COL},inplace=True)

    #merge original_df with na corrected df to get specific column information
    joined= pd.merge(df_long, original_df, on=[cons.NAME_COL, cons.FORMULA_COL, cons.LABEL_COL, cons.SAMPLE_COL],
                                                                                 how='left').fillna(0)     
    joined.loc[joined.Original_label == 0, cons.ORIGINAL_LABEL_COL] = joined.Label
    joined.drop(cons.LABEL_COL, axis=1, inplace=True)
    joined.rename(index= str, columns={cons.ORIGINAL_LABEL_COL :cons.LABEL_COL}, inplace=True)  
    return joined, eleme_corr_dict

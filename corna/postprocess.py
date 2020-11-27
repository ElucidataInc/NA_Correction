import pandas as pd
import numpy as np

from corna import constants as const
from corna import helpers as hlp
from isotopomer import Infopacket

def calculate_pool_total(df):
    """
    This function calculates pool total
    for all the metabolites in the input data file

    Args:
        df: dataframe for which pool total has to be calculated.
    Returns:
        final_df: dataframe which consists of the calculated values.
    """
    df[const.POOL_TOTAL_COL] = df[const.NA_CORRECTED_WITH_ZERO]
    df1 = df.groupby([const.SAMPLE_COL, const.NAME_COL])[const.POOL_TOTAL_COL].sum().reset_index()
    df.drop(const.POOL_TOTAL_COL, axis=1, inplace=True)
    df = df.merge(df1, on=[const.SAMPLE_COL, const.NAME_COL])
    return df


def pool_total(na_corr_df, colname):
    """
    This function calculates the pool total for each metabolite in a sample
    Args:
        na_corr_df: data frame with corrected intensities
        colname: Name of the column that contains corrected intensities
 
     Returns: grouped data frame with pool total
 
     """
    pool_total_df = na_corr_df.groupby(([NAME, SAMPLE]))
    pool_total_df = pool_total_df.apply(lambda x: x[x[colname] >= 0][colname].sum())
    return pool_total_df


def pool_total_MSMS(na_corr_df, colname):
    """
    This function calculates the pool total for each metabolite in a sample
    Args:
        na_corr_df: data frame with corrected intensities
        colname: Name of the column that contains corrected intensities

    Returns: grouped data frame with pool total
    """
    na_corr_df[METABOLITE_NAME] = na_corr_df[NAME].apply(get_metabolite)
    pool_total_df = na_corr_df.groupby(([METABOLITE_NAME, SAMPLE]))
    pool_total_df = pool_total_df.apply(lambda x: x[x[colname] >= 0][colname].sum())
    return pool_total_df

def fractional_enrichment(df):
    """
    This function calculates fractional enrichment
    for all the metabolites in the input data file
    Args:
        df: dataframe for which fractional enrichment has to be calculated.
    Returns:
        final_df: dataframe which consists of the calculated values.
    """
    final_df = pd.DataFrame()
    df = df.filter([const.SAMPLE_COL, const.NAME_COL, const.LABEL_COL, const.FORMULA_COL, const.NA_CORRECTED_COL])
    df = hlp.replace_negatives_in_column(df,const.NA_CORRECTED_WITH_ZERO, const.NA_CORRECTED_COL)

    df = calculate_pool_total(df)
    df[const.FRACTIONAL_ENRICH] = df[const.NA_CORRECTED_WITH_ZERO]/df[const.POOL_TOTAL_COL]
    final_df = df.fillna(0)
    if const.NA_CORRECTED_COL in final_df.columns:
        final_df.drop([const.NA_CORRECTED_COL], axis=1, inplace=True)
    if const.NA_CORRECTED_WITH_ZERO in final_df.columns:
        final_df.drop([const.NA_CORRECTED_WITH_ZERO], axis=1, inplace=True)
    return final_df


def zero_if_negative(num):
     """
     This function replaces negative numbers by zero_if_negative

     Args:
         num : any int value

     Return:
         1. zero if num is negative
         2. number itself if num is non negative
     """
     return 0 if num < 0 else num

def sum_intensities(fragments_dict):
    """
    This function calculates the sum of corrected intensities for a given sample

    Args:
        fragments_dict : dictionary of the form, example : {'Aceticacid_C13_1': [C2H4O2,
                         {'sample_1': array([ 0.0164])}, False, 'Aceticacid']
    Returns:
        sum_dict :  dictionary of sum of all corrected intensities for each sample
    """
    all_frag_info = fragments_dict.values()

    sample_names = []
    for frag in all_frag_info:
        sample_names.extend(frag.data.keys())
    sample_names = list(set(sample_names))

    sum_dict = {}

    for sample_name in sample_names:
        sum_dict[sample_name] = sum(value.data.get(sample_name,0) for value in all_frag_info)

    return sum_dict

def enrichment(fragments_dict, decimals):
    """
    This function calculates the fractional enrichment for each label
    Fractional enrichment[sample1] = Corrected Intensity/ Sum of corrected intensities of all labels
 
     Args:
         fragments_dict : dictionary of the form, example : {'Aceticacid_C13_1': [C2H4O2,
                          {'sample_1': array([ 0.0164])}, False, 'Aceticacid']
 
         decimals : number of significant digits to keep

     Returns:
         fragments_fractional : fragment dictionary model of fractional enrichment values
     """
    fragments_fractional = {}
    sum_dict = sum_intensities(fragments_dict)
 
    for key, value in fragments_dict.iteritems():
        fractional_data = {}
        for sample_name, intensity in value.data.iteritems():
            if not sum_dict[sample_name] == 0:
                fractional_data[sample_name] = np.around(
                     intensity / sum_dict[sample_name], decimals)
            else:
                fractional_data[sample_name] = 0
                warnings.warn("{} {} {} {}".format('sum of labels is zero for sample ', sample_name.encode('utf-8'),
                                                    ' of ', (value.name).encode('utf-8')))
        fragments_fractional[key] = Infopacket(
            value.frag, fractional_data, value.unlabeled, value.name)
 
    return fragments_fractional

def replace_vals(sample_int_dict):
     """
     This function replace negatives by zero in sample intensity dictionary
     Args:
         sample_int_dict : dictionary with keys as samples and values as corrected intensities
 
     Returns:
         dict_replaced_vals : sample_int_dict with negative intensities replaced by zeroes
     """
     dict_replaced_vals = {}
 
     for sample, intensity in sample_int_dict.iteritems():
        dict_replaced_vals[sample] = zero_if_negative(intensity)
 
     return dict_replaced_vals


def replace_negative_to_zero(corrected_dict):
     """
     This function replaces negative intensity values by zero from list of intensity
     in the standardised model dictionary
 
     Args:
         corrected_dict : nested dictionary (std model) with NA corrected intensity values
 
     Returns:
         post_proc_dict : returns nested dictionary with negative values replaced
 
     """
 
     post_proc_dict = {}
 
     for frag_key, frag_info in corrected_dict.iteritems():
         sample_int_dict = frag_info.data
         dict_replaced_vals = replace_vals(sample_int_dict)
         post_proc_dict[frag_key] = Infopacket(
             frag_info.frag, dict_replaced_vals, frag_info.unlabeled, frag_info.name)
 
     return post_proc_dict

def replace_negatives(na_corr_dict):
     """
     This function is a wrapper around replace_negatives_to_zero, it performs this
     function for all metabolites in the input data file
 
     Args:
         corrected_dict : nested dictionary (std model) with NA corrected intensity values
 
     Returns:
         post_proc_dict : returns nested dictionary with negative values replaced for
                          all the metabolites in the input data file
 
     """
     post_processed_dict = {}
     for metabolite, fragment_dict in na_corr_dict.iteritems():
         post_processed_dict[
             metabolite] = replace_negative_to_zero(fragment_dict)
 
     return post_processed_dict

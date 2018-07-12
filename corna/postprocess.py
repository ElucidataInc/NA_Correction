import pandas as pd

from corna import constants as const


def replace_negatives(df):
    """
    This function replaces negative intensity values to zero.
    """
    df[const.NA_CORRECTED_WITH_ZERO] = df[const.NA_CORRECTED_COL].clip(lower=0)
    return df


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
    df = replace_negatives(df)

    df = calculate_pool_total(df)
    df[const.FRACTIONAL_ENRICH] = df[const.NA_CORRECTED_WITH_ZERO]/df[const.POOL_TOTAL_COL]
    final_df = df.fillna(0)
    if const.NA_CORRECTED_COL in final_df.columns:
        final_df.drop([const.NA_CORRECTED_COL], axis=1, inplace=True)
    if const.NA_CORRECTED_WITH_ZERO in final_df.columns:
        final_df.drop([const.NA_CORRECTED_WITH_ZERO], axis=1, inplace=True)
    return final_df

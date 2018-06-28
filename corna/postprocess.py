import pandas as pd

def replace_negatives(df):
    """
    This function replaces negative intensity values to zero.
    """
    df['NA Corrected with zero']= df['NA Corrected'].clip(lower=0)
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
    df['Pool_total']=df['NA Corrected with zero']
    df1= df.groupby(['Sample','Name'])['Pool_total'].sum().reset_index()
    df.drop('Pool_total', axis=1, inplace=True)
    df= df.merge(df1, on=['Sample', 'Name'])
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
    final_df= pd.DataFrame()
    df= df.filter(['Sample', 'Name', 'Label','Formula', 'NA Corrected'])
    df= replace_negatives(df)

    df= calculate_pool_total(df)
    df['Fractional_enrichment']= df['NA Corrected with zero']/df['Pool_total']
    final_df= df.fillna(0)
    if 'NA Corrected' in final_df.columns:
        final_df.drop(['NA Corrected'], axis=1, inplace=True)
    if 'NA Corrected with zero' in final_df.columns:
        final_df.drop(['NA Corrected with zero'], axis=1, inplace=True)
    return final_df

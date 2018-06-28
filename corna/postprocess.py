import pandas as pd

def replace_negatives(df):
    """
    This function replaces negative intensity values to zero.
    """
    df['NA Corrected with zero']= df['NA Corrected'].clip(lower=0)
    return df

def pool_total(df):
    """
    This function calculates the pool total for each metabolite in a sample
    Args:
        df: data frame with corrected intensities
    Returns: df: datframe which column consisting of pool total for each metabolite
             total: pool total for a single metabolite.

    """
    total= df['NA Corrected with zero'].sum()
    df['Pool_Total']= total
    return df, total

def fractional_enrichment(df):
    """
    This function calculates fractional enrichment
    for all the metabolites in the input data file

    Args:
        df: dataframe for which fractional enrichment has to be calculated.
    Returns:
        final_df: dataframe which consists of the calculated values.
    """
    df= df.filter(['Sample', 'Name', 'Label','Formula', 'NA Corrected'])
    df= replace_negatives(df)
    final_df= pd.DataFrame()

    for samp in df.Sample.unique():
        samp_df= df[df['Sample']== samp]
        for metab in samp_df.Name.unique():
            metab_df= samp_df[samp_df['Name']== metab]
            metab_df, pooltotal= pool_total(metab_df)
            metab_df['Fractional_enrichment']= metab_df['NA Corrected with zero']/pooltotal
            final_df= final_df.append(metab_df)
    final_df= final_df.fillna(0)
    if 'NA Corrected' in final_df.columns:
        final_df.drop(['NA Corrected'], axis=1, inplace=True)
    if 'NA Corrected with zero' in final_df.columns:
        final_df.drop(['NA Corrected with zero'], axis=1, inplace=True)
    return final_df



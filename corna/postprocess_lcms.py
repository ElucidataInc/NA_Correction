"""
Module to calculate fractional enrichment if dataframe is in wide format.

"""

import pandas as pd

def replace_negatives(df):
    df['NA Corrected with zero']= df['NA Corrected'].clip(lower=0)
    return df

def pool_total(df, samp):
    total= df[samp].sum()
    return total

def fractional_enrichment_lcms(df):
    #df=pd.read_csv('/home/priyanka/Downloads/Test/final_out_fc.csv', index_col= 0)
    df= replace_negatives(df)
    original_df= df
    pool_dict={}
    metab_pool_dict={}
    df=df.pivot_table(index=['Name','Formula','Label'], columns='Sample', values='NA Corrected with zero')
    df = df.rename_axis(None, axis=1).reset_index()

    sample_df = df.select_dtypes(include=['float64', 'int'])
    sample_list= sample_df.columns.tolist()

    final_df= pd.DataFrame()
    output= pd.DataFrame()

    for metab in df.Name.unique():
        metab_df=df[df['Name']==metab]
        for sample in sample_list:
            total= pool_total(metab_df, sample)
            pool_dict[sample]= total
            metab_df[sample]= metab_df[sample]/total
        metab_pool_dict[metab]= pool_dict
        final_df=final_df.append(metab_df)

    df_long = pd.melt(final_df, id_vars=['Name', 'Formula', 'Label'])
    df_long.rename(columns={'variable': 'Sample', 'value':'Fractional Enrichment'},inplace=True)
    
    for metab in df_long.Name.unique():
        for sample in sample_list:
            df_long['Pool Total']= metab_pool_dict[metab][sample]
            output= output.append(df_long)
    final_df.to_csv('/home/priyanka/Downloads/Test/fractional_lcms1.csv')
    return final_df

df= pd.read_csv('/home/priyanka/Downloads/Test/lcms_nacorr.csv', index_col=0)
fractional_enrichment_lcms(df)
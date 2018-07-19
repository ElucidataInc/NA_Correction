import os
import pandas as pd
import pickle

import corna

path_dir = os.path.join(os.path.dirname(__file__),"test_data")
gcms_corr_merged = pickle.load(open(os.path.join(path_dir,"gcms_corr.p"),"r"))

def test_gcms_full():
    input_file = os.path.join(path_dir,"GCMS_raw.xlsx")
    maven_data = corna.read_file(input_file)
    iso_tracers = ['C13']
    eleme_corr = {'C': ['R','O','H']}
    na_dict = corna.get_na_value_dict()
    na_corr_df = corna.na_correction(maven_data, iso_tracers, 10 , na_dict, eleme_corr, autodetect=False)
    postprocessed_out = corna.replace_negatives(na_corr_df)
    frac_enrichment = corna.fractional_enrichment(postprocessed_out)
    df_list = [na_corr_df, frac_enrichment, postprocessed_out, maven_data]
    merged_results_df = corna.merge_multiple_dfs(df_list)
    print gcms_corr_merged
    assert pd.DataFrame.equals(merged_results_df, gcms_corr_merged)

test_gcms_full()
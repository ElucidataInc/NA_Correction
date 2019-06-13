
# coding: utf-8

# In[48]:

import pandas as pd
import numpy as np
import re
import json

import corna.constants as const
from corna.helpers import get_isotope_na, replace_negatives_in_column, merge_multiple_dfs
from corna.postprocess import fractional_enrichment
#from corna.algorithms.background_correction import background_correction
from corna.algorithms.background_correction import background_correction
from corna.inputs.multiquant_parser import mq_df_to_fragmentdict
from corna.algorithms.mimosa_nacorr import na_correction_mimosa
from corna.inputs import multiquant_parser
from corna.output import convert_to_df
from corna.postprocess import fractional_enrichment, replace_negatives
# In[49]:


print "in ==============="
raw_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/Raw_intensity_file.csv')
#print(raw_df.head())
print "\nmetadata mq"
metadata_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/metadata_mq_file.csv')
#print(metadata_df.head())
print "\n metadata sample file"
sample_metadata = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/sample_metadata.csv')
#print(sample_metadata.head())


# raw_df= pd.read_csv('raw_intensity_file.csv')
# print(raw_df.head())
# #metadata_df= pd.read_excel('metadata_mq.xlsx')
# metadata_df= pd.read_csv('metadata_mq_file.csv')
# print(metadata_df.head())
# #sample_metadata = pd.read_excel('metadata_sample.xlsx')
# sample_metadata = pd.read_csv('sample_metadata.csv')
# print(sample_metadata.head())
isBackground= False
isotope_dict= const.ISOTOPE_NA_MASS

# In[50]:

msms_df, list_of_replicates, sample_background = multiquant_parser.merge_mq_metadata(raw_df, metadata_df, sample_metadata)
print(msms_df.head())
msms_df.to_csv("msms_df.csv")
# print(list_of_replicates)
# print(sample_background)
# In[51]:

#background_corr_dict = multiquant_parser.mq_df_to_fragmentdict(msms_df)
#json = json.dumps(str(background_corr_dict))
#f = open("bg.json","w")
#f.write(json)

print("printing list of replacicates and sample background")
#print(list_of_replicates)
#print(sample_background)
#background_corr  = background_correction(background_corr_dict, list_of_replicates, sample_background)
background_corr  = background_correction(msms_df, list_of_replicates, sample_background)
bg_corrected_df = convert_to_df(background_corr, True, colname='Background Corrected')
bg_corrected_df.to_csv("background.csv")

na_corrected = na_correction_mimosa(bg_corrected_df, msms_df, const.ISOTOPE_NA_MASS)
na_corr_df = convert_to_df(na_corrected, True, colname='NA Corrected')
print(na_corr_df.head())
na_corr_df = replace_negatives_in_column(na_corr_df,const.NA_CORRECTED_WITH_ZERO, const.NA_CORRECTED_COL)
print(na_corr_df.head())
na_corr_df.to_csv('nacorrected.csv')
#print na_corrected

#merged_df = merge_multiple_dfs([bg_corrected_df, msms_df])
#merged_df.to_csv('merged_df_1.csv', index=False)
#metabolite_dict = mq_df_to_fragmentdict(merged_df, 'Background Corrected')
#print("metabolite dict")
#print(metabolite_dict)
#nacorrected = na_correction_mimosa(metabolite_dict)
#na_corr_df = convert_to_df(nacorrected, True, colname='NA corrected')
#na_corr_df.to_csv('nacorrected.csv')

# postprocessed_out = replace_negatives(nacorrected)
# postprocessed_out_df = convert_to_df(postprocessed_out, True, colname='NA_Corrected_with_zero')
# postprocessed_out_df.to_csv("postprocessed_out.csv")

# calculate fractional enrichment on post processed data
frac_enrichment = fractional_enrichment(na_corrected)
frac_enr_df = convert_to_df(frac_enrichment, True, colname='Frac Enrichment')
frac_enr_df = merge_multiple_dfs([msms_df, bg_corrected_df, na_corr_df, frac_enr_df])
frac_enr_df.to_csv("frac_enr.csv")

# 
# # In[53]:
# 
# fractional_enriched_df = fractional_enrichment(na_corrected)
# print fractional_enriched_df
# 
# # In[54]:
# 
# merged_df = merge_multiple_dfs([bg_corrected_df,na_corrected,fractional_enriched_df])
# print merged_df
# 
# # In[ ]:
# 


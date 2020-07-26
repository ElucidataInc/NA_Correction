
# coding: utf-8

# This notebook can be used to calculate NA Corrected intensities as well as fractional enrichment for LCMS input assuming high resolution. The example shows a dataset having two labels C13 and N15:
# 
#  - C13N15_lcms_high_res.csv - demo raw MS intensity file containing intensities for C17H27N3O17P2 from the data file of the repository published by Carreer William et al. in 2013

# In[1]:


import pandas as pd
import numpy as np
import re

from corna.inputs import maven_parser as parser
import corna.constants as const
from corna.helpers import replace_negatives_in_column, merge_multiple_dfs
from corna.algorithms.nacorr_lcms import na_correction
from corna.postprocess import fractional_enrichment


# Reading raw file and merging with sample metadata if present

# In[2]:


#raw_df = pd.read_csv('C13N15_lcms_high_res.csv')
#sample_metadata = pd.read_csv('meta_sample_lcms_high_res.csv')

raw_df = pd.read_csv('C13_lcms_high_res_newalgo.csv')
#if sample metadata not present, set it to empty dataframe
sample_metadata = pd.DataFrame()

merged_df, iso_tracer_data, element_list = parser.read_maven_file(raw_df, sample_metadata)
merged_df.head()


# Dictionary containing natural abundance values for the common isotopes found in nature. It can be defined by the user or one can use the default values from the package. The format of the dictionary is as shown below: 
# 
# {E:[M0, M1, ..Mn]} where E is the element symbol and the natural abundance fraction is in the increasing order of masses. For example:

# In[3]:


#user defined
# na_dict = {'O': [0.99757, 0.00038, 0.00205], 'H':[0.99985, 0.00015], 'N': [0.99632, 0.00368], 
#             'C': [0.9892, 0.0108], 'Si':[0.922297, 0.046832, 0.030872], 'S':[0.9493, 0.0076, 0.0429, 0, 0.0002]}


# Performing na_correction and inputs not relevant for this workflow are set as empty, using default dictionary from the package

# In[ ]:


na_corr_df, ele_corr_dict = na_correction(merged_df, iso_tracers=['C13'], res_type='ultra high res')

#for user defined NAdictionary
#na_corr_df, ele_corr_dict = na_correction(merged_df, iso_tracers=['C13', 'N15'], eleme_corr={}, na_dict=na_dict)

na_corr_df = replace_negatives_in_column(na_corr_df, const.NA_CORRECTED_WITH_ZERO, const.NA_CORRECTED_COL)
na_corr_df


# Calculating fraction enrichments, merging all data into a single file and saving as 'C13N15_lcms_high_res_corrected.csv'

# In[4]:


frac_enr_df = fractional_enrichment(na_corr_df)
frac_enr_df


# In[5]:


output_df = merge_multiple_dfs([merged_df, na_corr_df, frac_enr_df])
output_df


# In[6]:


merged_df.to_csv('C13N15_lcms_high_res_corrected.csv')


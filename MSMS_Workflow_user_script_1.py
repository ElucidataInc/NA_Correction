
# coding: utf-8

# In[48]:

import pandas as pd
import numpy as np
import re

import corna.constants as const
from corna.helpers import get_isotope_na, replace_negatives_in_column, merge_multiple_dfs
from corna.postprocess import fractional_enrichment
from corna.algorithms.background_correction import background_correction
from corna.algorithms.nacorr_msms import na_correction_mimosa
from corna.inputs import multiquant_parser

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


isBackground= False
isotope_dict= const.ISOTOPE_NA_MASS

# In[50]:

msms_df, list_of_replicates = multiquant_parser.merge_mq_metadata(raw_df, metadata_df, sample_metadata)

# In[51]:

bg_corrected_df = background_correction(msms_df, list_of_replicates, isotope_dict=const.ISOTOPE_NA_MASS)
bg_corrected_df.to_csv('background.csv')
bg_corrected_df = replace_negatives_in_column(bg_corrected_df,const.BACKGROUND_WITH_ZERO, const.BACKGROUND_CORRECTED)
#print bg_corrected_df

# In[52]:

na_corrected = na_correction_mimosa(bg_corrected_df, True, const.ISOTOPE_NA_MASS)
na_corrected = replace_negatives_in_column(na_corrected,const.NA_CORRECTED_WITH_ZERO, const.NA_CORRECTED_COL)
#print na_corrected

# In[53]:

fractional_enriched_df = fractional_enrichment(na_corrected)
#print fractional_enriched_df
fractional_enriched_df.to_csv('frac_enr_nacor.csv')

# In[54]:

merged_df = merge_multiple_dfs([bg_corrected_df,na_corrected,fractional_enriched_df])
#print merged_df

# In[ ]:



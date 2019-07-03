import os
import corna

# path to directory where multiquant text data files are present
# path_dir = '/Users/raa/OneDrive/Elucidata_Sini/NA_correction/Demo/data_yale/'
# path_dir = os.path.join(os.path.dirname(__file__), 'data_yale')
# 
# # read multiquant data files and combine them
# mq_files = corna.read_multiquant(os.path.join(path_dir,'tca_data_mq'))
# 
# # read multiquant metadata file
# mq_metadata = corna.read_multiquant_metadata(os.path.join(path_dir, 'mq_metadata.xlsx'))
# mq_sample_metadata = corna.read_multiquant_metadata(os.path.join(path_dir, 'metadata_samples.xlsx'))


from corna.inputs.multiquant_parser import merge_mq_metadata, mq_df_to_fragmentdict
from corna.algorithms.mimosa_bgcorr import met_background_correction
from corna.algorithms.mimosa_nacorr import na_correction_mimosa
from corna.postprocess import fractional_enrichment, replace_negatives
from corna.output import convert_to_df, convert_to_df_nacorr_MSMS
from corna.helpers import merge_multiple_dfs
import pandas as pd

print "in ==============="
sample_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/Raw_intensity_file.csv')
print(sample_df.head())
print "\nmetadata mq"
metadata_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/metadata_mq_file.csv')
print(metadata_df.head())
print "\n metadata sample file"
metadata_sample = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/sample_metadata.csv')
print(metadata_sample.head())

print "\n # merge multiquant files and metadata files"
merge_mq_metdata, list_of_replicates, sample_background = corna.merge_mq_metadata(sample_df, metadata_df, metadata_sample)
merge_mq_metdata.to_csv('merged_input_data.csv')

background_corr_dict = corna.mq_df_to_fragmentdict(merge_mq_metdata)
background_corr = corna.met_background_correction(background_corr_dict, list_of_replicates, sample_background)
# # convert background noise corrected dictionary to dataframe
background_corr_df = corna.convert_to_df(background_corr, True, colname='Background Corrected')
background_corr_df.to_csv("background.csv")

merged_df = corna.merge_multiple_dfs([background_corr_df, merge_mq_metdata])
metabolite_dict = corna.mq_df_to_fragmentdict(merged_df, 'Background Corrected')
nacorrected = corna.na_correction_mimosa(metabolite_dict)
na_corr_df = corna.convert_to_df(nacorrected, True, colname='NA corrected')
na_corr_df.to_csv('nacorrected.csv')
print(na_corr_df.head())

# # Replace negative values by zero on NA corrected data - optional
merged_df = corna.merge_multiple_dfs([na_corr_df, merge_mq_metdata])
metabolite_dict = corna.mq_df_to_fragmentdict(merged_df, 'NA corrected')
enrichment_dict = corna.na_correction_mimosa(metabolite_dict)

postprocessed_out = corna.replace_negatives(nacorrected)
postprocessed_out_df = corna.convert_to_df(postprocessed_out, True, colname='NA_Corrected_with_zero')
postprocessed_out_df.to_csv("postprocessed_out.csv")

# # calculate fractional enrichment on post processed data
frac_enrichment = corna.fractional_enrichment(postprocessed_out)
frac_enr_df = corna.convert_to_df(frac_enrichment, True, colname='Frac Enrichment')
frac_enr_df = corna.merge_multiple_dfs([merge_mq_metdata, background_corr_df, na_corr_df, postprocessed_out_df, frac_enr_df])
frac_enr_df.to_csv("frac_enr_cornagithub.csv")

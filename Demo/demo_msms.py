import os
import pandas as pd
import corna

# path to directory where multiquant text data files are present
#path_dir = '/Users/raa/OneDrive/Elucidata_Sini/NA_correction/Demo/data_yale/'
path_dir = os.path.join(os.path.dirname(__file__), 'data_msms')


# read multiquant data files and combine them
mq_files_fol = (os.path.join(path_dir,'tca_data_mq'))

# read multiquant metadata file
mq_raw_file = corna.read_file(os.path.join(mq_files_fol, 'raw_intensity_file.csv'))
mq_metadata = corna.read_file(os.path.join(mq_files_fol, 'metadata.xlsx'))

# get sample metadata file for background Correction
mq_sample_metadata = corna.read_file(os.path.join(mq_files_fol, 'metadata_sample.xlsx'))

# merge multiquant files and metadata files
merge_mq_metdata, list_of_replicates = corna.merge_mq_metadata(mq_raw_file, mq_metadata, mq_sample_metadata)

#merge_mq_metdata.to_csv(mq_files_fol + 'merged_input_data.csv')

# perform NA correction
nacorrected = corna.na_correction_mimosa(merge_mq_metdata,False)
#print nacorrected

postprocessed_out = corna.replace_negatives(nacorrected)
# print postprocessed_out

# # calculate fractional enrichment on post processed data
frac_enrichment = corna.fractional_enrichment(postprocessed_out)
frac_enr_df = corna.merge_multiple_dfs([merge_mq_metdata, nacorrected, frac_enrichment])
#print frac_enr_df
#
# # save any dataframe at given path
#frac_enr_df.to_csv(mq_files_fol + 'frac_enrichment_1.csv')

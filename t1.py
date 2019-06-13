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


#from corna.inputs.multiquant_parser import merge_mq_metadata, mq_df_to_fragmentdict
#from corna.algorithms.mimosa_bgcorr import met_background_correction
#from corna.algorithms.mimosa_nacorr import na_correction_mimosa
#from corna.postprocess import fractional_enrichment, replace_negatives
#from corna.output import convert_to_df, convert_to_df_nacorr_MSMS
#from corna.helpers import merge_multiple_dfs
import pandas as pd

print "in ==============="
sample_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/Raw_intensity_file.csv', header=0)
print(sample_df.head())

sample_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/Raw_intensity_file.csv', header=head)
print(sample_df.head())

# print "\nmetadata mq"
# metadata_df = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/metadata_mq_file.csv')
# print(metadata_df.head())
# print "\n metadata sample file"
# metadata_sample = pd.read_csv('/home/elucidata/work_station/isocorrect_shef/files/sample_metadata.csv')
# print(metadata_sample.head())


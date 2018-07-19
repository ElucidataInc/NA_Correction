import os

import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest

from corna.inputs import multiquant_parser

DIR_PATH = os.path.dirname(os.path.realpath(__file__))
MQ_FILE_PATH = os.path.join(DIR_PATH, 'test_input_validation_data', "raw_mq.txt")
MQ_METADATA_PATH = os.path.join(DIR_PATH, 'test_input_validation_data', "metadata_mq.xlsx")
MQ_SAMPLE_METADATA_PATH = os.path.join \
    (DIR_PATH, 'test_input_validation_data', "metadata_sample.xlsx")
mq_df= df = pd.read_csv(MQ_FILE_PATH, sep='\t|,', engine='python')
metadata_df= pd.read_excel(MQ_METADATA_PATH)
sample_metadata_df= pd.read_excel(MQ_SAMPLE_METADATA_PATH)

INPUT_FILES = {"mq_file_path": MQ_FILE_PATH, \
                   "mq_metadata_path": MQ_METADATA_PATH, \
                   "mq_sample_metadata_path": MQ_SAMPLE_METADATA_PATH \
                   }
INPUT_FILES_WITHOUT_METADATA = {"mq_file_path": MQ_FILE_PATH,
                                "mq_metadata_path": MQ_METADATA_PATH,
                                "mq_sample_metadata_path": ''}

raw_df = pd.DataFrame({'Label': { 0: 'O18_191.0_111.0', 1: 'O18_193.0_113.0', 2: 'O18_195.0_115.0'},
                        'Parent Formula':{ 0: 'C6H7O7', 1: 'C6H7O7', 2: 'C6H7O7'},
                        'Formula': {0: 'C5H3O3',1: 'C5H3O3', 2: 'C5H3O3'}} )

mass_df = pd.DataFrame({
    'iso_mass': {0: 193.0, 1: 191.0},
    'mol_mass': {0: 191, 1:191}
})
 

def test_add_mass_and_no_of_atoms_info_frm_label():
    result_df, isotracer = multiquant_parser.add_mass_and_no_of_atoms_info_frm_label(raw_df)
    print result_df
    No_of_labeled_atoms_parent = [0, 1, 2]
    No_of_labeled_atoms_daughter = [0, 1, 2]
    assert result_df['Number of labeled atoms parent'].tolist() == No_of_labeled_atoms_parent
    assert result_df['Number of labeled atoms daughter'].tolist() == No_of_labeled_atoms_daughter

def test_get_num_labeled_atoms():
    result = multiquant_parser.get_num_labeled_atoms('O18', mass_df['iso_mass'], mass_df['mol_mass'])
    label_no= [1, 0]
    assert result.tolist() == label_no

def test_get_parent_daughter_iso_mass_col_and_isotracer_from_label():
    result_df = multiquant_parser.get_parent_daughter_iso_mass_col_and_isotracer_from_label(raw_df)
    isotracer = ['O18', 'O18', 'O18']
    iso_mass_parent = [191.0, 193.0, 195.0]
    assert result_df['Isotopic Tracer'].tolist() == isotracer
    assert result_df['isotopic_mass_parent'].tolist() == iso_mass_parent

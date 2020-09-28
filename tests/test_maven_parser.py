import os
import pytest
import pandas as pd

from pandas.util.testing import assert_frame_equal

from corna import custom_exception
from corna.inputs import maven_parser
import constants
from fixtures import *


df= pd.DataFrame({'Name': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic' , 3: 'Acetic'},
                        'Intensity': {0: 0.3624, 1: 0.040349999999999997, 2: 0.59724999999999995, 3: 0.11},
					    'Formula': {0: 'H4C2O2', 1: 'H4C2O2', 2: 'H4C2O2', 3: 'H4C2'},
					    'info2': {0: 'culture_1', 1: 'culture_1', 2: 'culture_1', 3: 'culure_1'},
					    'Sample': {0: 'sample_1', 1: 'sample_1', 2: 'sample_1', 3: 'sample_1'},
                        'Label': {0: 'C13N15-label-1-1', 1: 'C12-PARENT', 
                                  2: 'C13N15-label-0-1', 3: 'N15-label-2'}

})


df_label= pd.DataFrame({'Label': {0: 'C13N15-label-1-1', 1: 'C12-PARENT', 
                                  2: 'C13N15-label-0-1', 3: 'N15-label-2'}

})

df_label_new= pd.DataFrame({'Label': {0: 'C13N15-label-1-1', 1: 'C13N15-label-0-0', 
                                  2: 'C13N15-label-0-1', 3: 'C13N15-label-0-2'}

})


df_isotracers= pd.DataFrame({'C13': {0: 1, 1: 0, 2: 0, 3: 0},
                              'N15': {0: 1, 1: 0, 2: 1, 3: 2}
                            })
iso_tracers= ['C13', 'N15']

def test_get_isotope_columns_frm_label_col():
    df = maven_parser.get_isotope_columns_frm_label_col(df_label, iso_tracers)
    assert df.equals(df_isotracers)


def test_check_duplicates_in_list():
    check_list= [1, 3, 4, 1]
    assert maven_parser.check_duplicates_in_list(check_list) == [1]

def test_check_duplicates_in_list_no_duplicate():
    check_list= [1, 3, 4]
    assert maven_parser.check_duplicates_in_list(check_list) == []

def test_check_duplicates_in_list_empty():
    check_list= []
    assert maven_parser.check_duplicates_in_list(check_list) == []

def test_create_label_col_frm_isotope_columns():
    df = maven_parser.create_label_column_frm_isotope_columns(df_isotracers, iso_tracers)
    assert df['Label'].equals(df_label_new['Label'])

def test_save_original_label_and_processed_label():
    result = maven_parser.save_original_label_and_processed_label(df, iso_tracers)
    assert result['Original_label'].equals(df['Label']) and result['Label'].equals(df_label_new['Label'])

def test_check_df_empty():
    assert maven_parser.check_df_empty(pd.DataFrame())


def test_check_error_present():
    logs = {constants.VALIDATION_ERROR: ['There is one erroe'],
            constants.VALIDATION_WARNING: []}
    assert maven_parser.check_error_present(logs)


def test_get_extracted_isotracer():
    assert maven_parser.get_extracted_isotracer('C13-label-1') == 'C13'
    assert maven_parser.get_extracted_isotracer('C12 PARENT') == 'C12 PARENT'


def test_get_extraced_isotracer_df(get_maven_df):
    maven_df = get_maven_df
    test_assert = ['C12 PARENT', 'C13', 'C13']
    assert list(maven_parser.get_extraced_isotracer_df(maven_df)) == test_assert


def test_isotracer_dict(get_maven_df):
    maven_df = get_maven_df
    assert maven_parser.get_isotracer_dict(maven_df) == {'C13': 2, 'C12 PARENT': 1}


def test_get_extracted_element():

    assert maven_parser.get_extracted_element('C3H2O6') == {'H': 2, 'C': 3, 'O': 6}
    assert maven_parser.get_extracted_element('SiH2O6') == {'H': 2, 'Si': 1, 'O': 6}
    assert maven_parser.get_extracted_element('C3H2KFe') == {'H': 2, 'C': 3, 'K': 1, 'Fe': 1}


def test_get_element_list():
    input_df = read_csv(constants.MAVEN_FILE)
    assert maven_parser.get_element_list(input_df) == ['C', 'H', 'O', 'N']




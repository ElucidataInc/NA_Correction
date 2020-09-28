import os

import pandas as pd
import pytest

import corna.helpers as help



def test_parse_formula():
    assert help.parse_formula('C2H5O2') == {'C':2, 'H':5, 'O':2}


def test_parse_formula_single_atom():
    assert help.parse_formula('CH2O') == {'C':1, 'H':2, 'O':1}  


def test_get_global_isotope_dict():
    assert help.get_global_isotope_dict()['amu']['C12'] == 12


def test_get_atomic_weight():
    assert help.get_atomic_weight('Al') == 26.981538


def test_get_atomic_weight_wildcard():
    with pytest.raises(KeyError) as err:
        help.get_atomic_weight('PP')
    assert err.value.message == 'Element doesnt exist'


def test_check_if_isotope_in_dict():
    assert help.check_if_isotope_in_dict('C13') == True


def test_check_if_isotope_in_dict_false():
    assert help.check_if_isotope_in_dict('C14') == False


def test_get_isotope_details():
    assert help.get_isotope_element('C13') == 'C'


def test_get_isotope_na():
    assert help.get_isotope_na('C13') == 0.0111


def test_get_isotope_na_keyerror():
    with pytest.raises(KeyError):
        help.get_isotope_na('PP')


def test_get_isotope_mass():
    assert help.get_isotope_mass('C13') == 13


def test_get_isotope_mass_keyerror():
    with pytest.raises(KeyError):
        help.get_isotope_mass('C23')

def test_mol_weight():
    assert help.get_mol_weight('C2H5') == 29.0611


def test_get_isotope_natural():
    assert help.get_isotope_natural('C13') == 'C12'


def test_get_isotope_natural_keyerror():
    with pytest.raises(KeyError):
        help.get_isotope_natural('PP')


def test_read_file():
    path = 'incorrectpath.xlsx'
    with pytest.raises(IOError):
        file_in = help.read_file(path)


def test_read_file_ext():
    path = 'no_or_invalid_extension'
    with pytest.raises(IOError):
        file_in = help.read_file(path)


def test_read_file_csv():
    input_file_path = os.path.join(os.path.dirname(__file__), 'test_data/Data_lcms.csv')
    assert list(help.read_file(input_file_path)) == ['Label', 'Name', 'Formula', 'Area']


def test_read_file_txt():
    input_file_path = os.path.join(os.path.dirname(__file__), 'test_data/test_data.txt')
    assert list(help.read_file(input_file_path)) == ['Original Filename', 'Sample Name', 'Glucose Concentration', 'Time',
                                                     'Component Name', 'Mass Info', 'Area', 'Modified']

def test_filter_df_error():
    df = pd.DataFrame({'col1': [1, 2, 3], 'col2': [2, 2, 2]})
    filter_dict = {'col_1': [1]}
    with pytest.raises(KeyError):
        filter_df = help.filter_df(df, filter_dict)


def test_get_na_value_dict_O():
    assert help.get_na_value_dict()['O'] == [0.9976, 0.0004, 0.002]


def test_filter_df():
    df = pd.DataFrame({'Name': ['a', 'b', 'c'], 'Formula': ['a', 'b', 'c']})
    colval = {'Name': ['a']}
    assert help.filter_df(df, colval)['Name'].count() == 1


def test_get_unique_val():
    df = pd.DataFrame({'Name': ['a', 'b', 'c'], 'Formula': ['a', 'b', 'c']})
    with pytest.raises(KeyError) as err:
        help.get_unique_values(df, 'NAME')
        assert err.value.message == 'Column NAME not found in dataframe'


def test_concatenate_dataframes_by_col():
    df_1= pd.DataFrame({'Name':['a', 'b'], 'Formula':['d', 'e']})
    df_2 = pd.DataFrame({'mode':['+', '+']})
    assert list(help.concatenate_dataframes_by_col([df_1, df_2])) == ['Formula', 'Name', 'mode']


def test_parse_polyatom():
    assert help.parse_polyatom('Si32') == ('Si', 32)


def test_merge_multiple_df():
    df_1= pd.DataFrame({'Name': ['a', 'b'], 'Formula': ['d', 'e'],
                        'Label': ['f','g'], 'Sample': [1,2]})
    df_2 = pd.DataFrame({'Name': ['a'], 'Label': ['f'], 'Sample': [1], 'Formula': ['d']})
    result_df = help.merge_multiple_dfs([df_1, df_2])
    assert result_df.equals(df_2)


def test_get_metabolite_fragment():
    assert help.get_metabolite('Citrate 191/111') == 'Citrate'


def test_first_sub_second():
    assert help.first_sub_second([1, 2], [3, 4]) == [1, 2]


def test_get_isotope_keyerror():
    assert not help.check_if_isotope_in_dict('Ind5')
    

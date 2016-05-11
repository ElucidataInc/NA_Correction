import pytest
import pandas as pd
import corna.helpers as help

def test_get_atomic_weight():
    assert help.get_atomic_weight('Al') == 26.981538

def test_get_atomic_weight_wildcard():
    with pytest.raises(KeyError) as err:
        help.get_atomic_weight('R')
    assert err.value.message == 'Element doesnt exist'

def test_read_file():
	path = 'incorrectpath.xlsx'
	with pytest.raises(IOError):
		file_in = help.read_file(path)

def test_read_file_ext():
	path = 'no_or_invalid_extension'
	with pytest.raises(IOError):
		file_in = help.read_file(path)



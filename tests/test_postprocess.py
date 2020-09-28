import pytest
import numpy as np
import pandas as pd

import corna.postprocess as postprocess
import constants


df = pd.DataFrame({'Name': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'},
					   'Parent': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'},
					   'Label': {0: 'C12 PARENT', 1: 'C13-label-1', 2: 'C13-label-2'},
					   'NA Corrected': {0: 1.23, 1: -2.0, 2: 1.77},
					   'Formula': {0: 'H4C2O2', 1: 'H4C2O2', 2: 'H4C2O2'},
					   'info2': {0: 'culture_1', 1: 'culture_1', 2: 'culture_1'},
					   'Sample': {0: 'sample_1', 1: 'sample_1', 2: 'sample_1'}})

df_with_zero= pd.DataFrame({'Name': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'},
					   'Label': {0: 'C12 PARENT', 1: 'C13-label-1', 2: 'C13-label-2'},
					   'NA Corrected': {0: 1.23, 1: -2.0, 2: 1.77},
					   'NA Corrected with zero': {0: 1.23, 1:0, 2:1.77},
					   'Formula': {0: 'H4C2O2', 1: 'H4C2O2', 2: 'H4C2O2'},
					   'Sample': {0: 'sample_1', 1: 'sample_1', 2: 'sample_1'}})

df_negative = pd.DataFrame({'NA Corrected': {0: 123, 1: -78.0}})
df_positive = pd.DataFrame({'NA Corrected': {0: 123, 1: 78.0}})


def test_replace_negatives_zero():
	result = postprocess.replace_negatives(df_negative)
	output_list = [123, 0]
	assert result['NA Corrected with zero'].tolist()== output_list
	

def test_replace_negatives():
	result = postprocess.replace_negatives(df_positive)
	output_list = [123, 78]
	assert result['NA Corrected with zero'].tolist()== output_list


def test_calculate_pool_total():
	result= postprocess.calculate_pool_total(df_with_zero)
	output_list=[3.0, 3.0, 3.0]
	assert result['Pool_total'].tolist()== output_list

def test_fractional_enrichment():
	result= postprocess.fractional_enrichment(df)
	output_list=[0.41, 0.00, 0.59]
	assert result['Fractional enrichment'].tolist()== output_list








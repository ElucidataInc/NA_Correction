import os
import pandas as pd
import pytest

from corna.algorithms.nacorr_lcms import na_correction

single_tracers = ['C13']
multi_tracers = ['C13', 'N15']
na_dict = {'H':[0.98,0.01,0.01], 'C': [0.95, 0.05], 'S': [0.922297, 0.046832, 0.030872], 'O':[0.95,0.03,0.02], 'N': [0.8, 0.2]}
na_dict_auto= {'C':[0.9889,0.0111],
'N':[0.9964,0.0036],
'O':[0.9976,0.0004,0.002],
'H':[0.99985,0.00015],
'S':[0.95,0.0076,0.0424],
}

def test_na_corr_single_tracer():
	df = pd.DataFrame({'Name': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'},
					   'Parent': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'},
					   'Label': {0: 'C12 PARENT', 1: 'C13-label-1', 2: 'C13-label-2'},
					   'Intensity': {0: 0.3624, 1: 0.040349999999999997, 2: 0.59724999999999995},
					   'Formula': {0: 'H4C2O2', 1: 'H4C2O2', 2: 'H4C2O2'},
					   'info2': {0: 'culture_1', 1: 'culture_1', 2: 'culture_1'},
					   'Sample': {0: 'sample_1', 1: 'sample_1', 2: 'sample_1'}})

	eleme_corr = {}

	na_corr_df, corr_dct = na_correction(df, ['C13'], '', na_dict, eleme_corr,
									       autodetect=False)
	print na_corr_df

	output_list = [0.4015512465373961, 0.0023185595567865968, 0.59613019390581734]

	assert na_corr_df['NA Corrected'].tolist() == output_list


def test_na_corr_single_trac_indist():

	df = pd.DataFrame({'Name': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'}, \
	 'Parent': {0: 'Acetic', 1: 'Acetic', 2: 'Acetic'}, \
	  'Label': {0: 'C12 PARENT', 1: 'C13-label-1', 2: 'C13-label-2'}, \
	   'Intensity': {0: 0.2274, 1: 0.4361, 2: 0.25405}, \
	    'Formula': {0: 'H4C2O2', 1: 'H4C2O2', 2: 'H4C2O2'}, \
	    'info2': {0: 'culture_1', 1: 'culture_1', 2: 'culture_1'}, \
	     'Sample': {0: 'sample_1', 1: 'sample_1', 2: 'sample_1'}})

	eleme_corr = {'C': ['H', 'O']}

	na_corr_df, corr_dict = na_correction(df, ['C13'], '', na_dict, eleme_corr, 
											autodetect=False)

	output_list = [0.3035536244690365, 0.48572975053068485, 0.1961957568543734]

	assert na_corr_df['NA Corrected'].tolist() == output_list

def test_na_corr_multi_trac():
	df = pd.DataFrame({'Name': {0: 'L-Methionine', 1: 'L-Methionine'},
					   'Label': {0: 'C12 PARENT', 1: 'C13-label-1'},
					   'Intensity': {0: 0.203405, 1: 0.050069999999999996},
					   'Formula': {0: 'C5H10NO2S', 1: 'C5H10NO2S'},
					   'Sample': {0: 'sample_1', 1: 'sample_1'}})
	eleme_corr = {}
	na_corr_df, corr_dict = na_correction(df, ['C13', 'N15'], '',
											na_dict, eleme_corr,
											autodetect=False)

	output_list = [0.3285894465447474, -0.06571788930894949, -0.005306330643483816, 
				0.0010612661286967635, -0.007153470034923123, 0.0014306940069846248, 
				0.0007418786567104569, -0.00014837573134209142, -2.8152095292076256e-05, 
				5.630419058415252e-06, 3.775722416858374e-07, -7.55144483371675e-08]
	assert na_corr_df['NA Corrected'].tolist() == output_list


def test_na_corr_multi_trac_indist():
	df = pd.DataFrame({'Name': {0: 'L-Methionine', 1: 'L-Methionine'},
					   'Label': {0: 'C12 PARENT', 1: 'C13-label-1'},
					   'Intensity': {0: 0.203405, 1: 0.050069999999999996},
					   'Formula': {0: 'C5H10NO2S', 1: 'C5H10NO2S'},
					   'Sample': {0: 'sample_1', 1: 'sample_1'}})
	eleme_corr = {'C': ['H']}
	na_corr_df, corr_dict =  na_correction(df, ['C13', 'N15'],'',
														   na_dict, eleme_corr,
														   autodetect=False)
	output_list = [0.40215440095785504, -0.08043088019157102, -0.04547876075722601,
				 0.009095752151445204, -0.04503063995239236, 0.009006127990478473,
				 0.002629393437085129, -0.000525878687417026, 0.002627198654786159,
				 -0.000525439730957232, 3.300655199399871e-05, -6.6013103987997436e-06]
	assert na_corr_df['NA Corrected'].tolist() == output_list

def test_ppm_nacorrection():
	df = pd.DataFrame({'Name': {0: 'Pyruvic acid', 1: 'Pyruvic acid'},
					   'Label': {0: 'C12 PARENT', 1: 'C13-label-1'},
					   'Intensity': {0: 287515.1, 1: 7354.636},
					   'Formula': {0: 'C3H4O3', 1: 'C3H4O3'},
					   'Sample': {0: 'sample_1', 1: 'sample_1'}})
	eleme_corr = {}
	na_corr_df, corr_dict = na_correction(df, ['C13'], 50, na_dict_auto, eleme_corr, 
											autodetect=True)

	output_list = [297484.35034520662, -2557.5824049265775, -57.580996614900414, 0.54905682625541885]
	assert list(na_corr_df['NA Corrected']) == output_list
	assert corr_dict == {'Pyruvic acid': {'C': ['H', 'O17']}}



"""
def test_autodetect_nacorr():
	df = pd.DataFrame({'Name': {0: 'Pyruvic acid', 1: 'Pyruvic acid'},
					   'Label': {0: 'C12 PARENT', 1: 'C13N15-label-3-1'},
					   'Intensity': {0: 0.3303, 1: 0.5065},
					   'Formula': {0: 'C10H17N3O6S', 1: 'C10H17N3O6S'},
					   'Sample': {0: 'sample_1', 1: 'sample_1'}})
	#eleme_corr = {}
	eleme_corr = {'C':['H','O','O'], 'N':['S']}
	na_corr_df, corr_dict = na_correction(df, ['C13','N15'], 10, na_dict_auto, eleme_corr, 
											autodetect=False)
	print na_corr_df
"""
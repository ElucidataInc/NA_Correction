"""This module helps to do validation using validation package olmonk"""

from corna.olmonk import basic_validation, data_validation

def get_validation_df(path, required_columns=None):
	"""takes path of the file, validates it and returns result

	Using path of the file, it instantiates the validation class
	and use its methods to validate file and returns result of
	that. If file not passes the validation checks, then it will
	raise an error, otherwise it will returns validated df.
	Args:
		path: file path for which validation is needed

	Returns: instance of validated df
	"""
	try:
		basic_validator = get_class_inst(basic_validation.BasicValidator, path, required_columns)
		basic_validation_result(basic_validator)
		return basic_validator
	except Exception as e:
		raise

def get_class_inst(validator_class, file_path, required_columns):
	"""Instantiates class with its argument and returns it"""

	return validator_class(file_path, required_columns)

def basic_validation_result(basic_validator):
	"""Takes instance of class, do basic validation and raise error
	if any check fails
	"""

	try:
		basic_validator.check_path_exist()
		basic_validator.check_file_empty()
		basic_validator.check_if_convert_to_df()
	except Exception as e:
		raise

def data_validation_raw_df(df):
	"""do datavalidtaion for raw_file_df and returns report_df

	It takes df of raw_mq file, creates an instance of DataValidation
	using this df. It then does validation related to file and returns
	the report_df.

	Args:
		df: raw_mq_file df

	Returns:
		report_df contains report of error & warning in this df
	"""
	# :TODO: update doc when this function will be updated
	try:
		raw_df_validator = data_validation.DataValidator(df)
		raw_df_validator.missing_data()
		raw_df_validator.numerical(['Area'])
		raw_df_validator.pattern_match('Mass Info', '\d+.0 \/ \d+.0')
		return raw_df_validator.get_report_df()
	except Exception as e:
		raise

def data_validation_metadata_df(df):
	"""do datavalidtaion for metadata_mq_df and returns report_df

	It takes df of metadata_mq file, creates an instance of DataValidation
	using this df. It then does validation related to file and returns
	the report_df.

	Args:
		df: metadata_mq_file df

	Returns:
		report_df contains report of error & warning in this df
	"""
	# :TODO: update doc when this function will be updated
	try:
		metadata_df_validator = data_validation.DataValidator(df)
		metadata_df_validator.missing_data()
		metadata_df_validator.chemical_formula(['Formula', 'Parent Formula'])
		metadata_df_validator.value_in_constant('Isotopic Tracer', ['C13', 'N15'])
		return metadata_df_validator.get_report_df()
	except Exception as e:
		raise




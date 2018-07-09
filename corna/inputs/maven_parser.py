import pandas as pd

from corna.helpers import get_isotope_element, first_sub_second, parse_formula, chemformula_schema, read_file, check_column_headers 
from column_conventions import maven as maven_constants
import corna.constants as cons

def convert_labels_to_std(df, iso_tracers):
    """
    This function breaks the labels C13N15-label-1-1 to form columns as C13 and N15 and 
    1 and 1 as values respectively for all labels in a column.
    """
    def process_label(label):
        """
        This function converts the labels C13N15-label-1-1 in the form
        C13_1_N15_1
        """
        if 'PARENT' in label:
            return '_'.join('{}_0'.format(t) for t in iso_tracers)
        else:
            formula, enums = label.split('-label-')
            isotopes = set(''.join(map(str, i))
            for i in chemformula_schema.parseString(formula))
            msg = """iso_tracers must have all isotopes from input data
                    Got: {!r}
                    Expected: {!r}
                  """.format(iso_tracers, isotopes.union(iso_tracers))
            assert set(isotopes).issubset(set(iso_tracers)), msg
            # The final label must have all iso_tracers
            # Use zeroes as default, else the number from given label
            inmap = {i: 0 for i in iso_tracers}
            inmap.update({i: n for i, n in zip(isotopes, enums.split('-'))})
            # The order is important, so we don't map on inmap directly
            return '_'.join("{}_{}".format(i, inmap[i]) for i in iso_tracers)
    
    df[cons.LABEL_COL] = [process_label(l) for l in df[cons.LABEL_COL]]
    
    s=df[cons.LABEL_COL].apply(lambda x: x.split('_'))
    i=0
    for iso in iso_tracers:
        df[iso]= s.apply(lambda x: x[i+1])
        df[iso]= df[iso].astype(int)
        i=i+2

    del df[cons.LABEL_COL]
    return df

def process_corrected_df_for_metab(df, metabolite, isotracers, required_column):
    """
    This function returns formula dict for the formula of metabolite.
    formula- C5H4N3
    formula_dict- {'C':5,'H':4, 'N':3}
    """
    metabolite_df=df[df[cons.NAME_COL]==metabolite]
    formula = metabolite_df.Formula.unique()
    formula_dict=parse_formula(formula[0])
    metabolite_df.set_index(isotracers, inplace=True) 
    required_df= metabolite_df.filter(required_column)
    return required_df, formula, formula_dict

def check_duplicates_in_list(given_list):
    """
    This function checks for duplicates in the given list.
    It iterates over the list and whenever a new element is found,
    we add that to first_occurrence

    :param given_list:
    :return: duplicate_list : list of all the duplicates in given list
    """
    first_occurrence = set()
    duplicate_list = set()
    first_occurrence_add = first_occurrence.add
    duplicate_list_add = duplicate_list.add
    for item in given_list:
        if item in first_occurrence:
            duplicate_list_add(item)
        else:
            first_occurrence_add(item)
    return list(duplicate_list)


def create_label_column(df, isotracers):
    """
    This function creates the label column back by joining information from multiple splited columns.
    Example:
            C13 N15 Label
            1   1   C13N15-Label-1-1 
    """
    df =df.rename_axis(None, axis=1).reset_index()
    isogroup=""
    df[cons.LABEL_COL]=""
    for i in isotracers:
        isogroup += i
        df[cons.LABEL_COL] +='-'+ df[i].astype(str)
        del df[i]
    df[cons.LABEL_COL]= isogroup +'-label'+ df[cons.LABEL_COL]

    s=df[cons.LABEL_COL].apply(lambda x: x.split('-'))

    return df

def add_info_to_final_df(info_df, metab, formula, iso_tracers):
    """
    Adds required columns back to the na corrected Dataframe 
    Required columns include : Label, Name, Formula
    """
    info_df=create_label_column(info_df, iso_tracers)
    info_df[cons.NAME_COL]= metab
    info_df[cons.FORMULA_COL] = formula
    return info_df

def save_original_df(df, iso_tracers):
    """
    This function modifies Label column from C13-Label-1 to C13N15-Label-1-0.
    And then returns back the original dataframe.
    """
    df= df.filter([cons.NAME_COL, cons.FORMULA_COL, cons.LABEL_COL, cons.SAMPLE_COL, cons.INTENSITY_COL])
    df[cons.ORIGINAL_LABEL_COL]=df[cons.LABEL_COL]
    df= convert_labels_to_std(df, iso_tracers)
    df= create_label_column(df, iso_tracers)
    df.drop('index', axis=1, inplace= True)
    return df

def check_error_present(logs):
    """
    This function checks if any error is present in the validation
    logs.
    :param logs: validation logs after all the validation checks
    :return: Boolean True if error is present
    """
    if logs[cons.VALIDATION_ERROR]:
        return True
    else:
        return False

def get_extracted_isotracer(label):
    """
    This function takes label as an argumnet and returns the
    iso-tracer present in label. This helps in counting iso-tracer
    present in label column value.
    for ex: label = 'C13N15-label-1-2
            returns = 'C13N15
    :param label: label column value
    :return: extracted iso-tracer value
    """
    if label == cons.UNLABELLED_LABEL:
        return cons.UNLABELLED_LABEL
    else:
        return label.split('-label-')[0]


def get_extracted_element(formula):
    """
    This function parse chemical formula and returns different elements
     in the formula.
    :param formula: Formula to be parsed
    :return: dict with keys containing different elements in the formula
    """
    # TODO: Similar to get formula can be modified for further advancement
    return parse_formula(formula)


def get_extraced_isotracer_df(maven_df):
    """
    This function extract iso-tracer information of
    all the values of label column.
    :param maven_df: df whose iso-tracer needs to get extracted
    :return: extracted iso-tracers for each label column value
    """
    return maven_df[maven_constants.LABEL].apply(get_extracted_isotracer)


def get_isotracer_dict(maven_df):
    """
    This function counts the iso-tracer in label column of
     the maven df. Dict with isotracer as key and its number as value
     is returned.
    :param maven_df: maven_df whose iso-tracer needs to be counted
    :return: dict of iso-tracer and its number
    """
    isotracer_df = get_extraced_isotracer_df(maven_df)
    return isotracer_df.value_counts().to_dict()

def convert_inputdata_to_stdfrom(input_df):
    """
    This function convert the input data file(maven format) into standard data model. It gives
    the same format as in merged data (maven + metadata). This function can be used if the user
    does not wish to input metadata file and proceed as it is with the input data file

    Args:
        input_df : MAVEN input data in the form of pandas dataframe
    Returns:
        std_form_df : dataframe with input data in standard data model format that can be
                    used in further processing
    """
    
    long_form = melt_df(input_df)
    std_form_df = column_manipulation(long_form)
    return std_form_df

def column_manipulation(df):
    """
    This function adds a parent column to dataframe and renames
    the sample and intensity columns in order to bring it in standard data model

    Args:
        df : dataframe on which column changes are to be done

    Returns:
        df : df with added columns and standard column names
    """
    df[cons.PARENT_COL] = df[maven_constants.NAME]
    df.rename(
        columns={cons.VAR_COL: cons.SAMPLE_COL, cons.VAL_COL: cons.INTENSITY_COL},
        inplace=True)

    return df

def check_df_empty(df):
    return df.empty

def get_element_list(maven_df):
    """
    This function gives list of uniqu element present in the formula
    column of maven df.
    :param maven_df:  df whose element is to be listed out
    :return: list of unique element in formula column
    """
    element_dict = dict()
    extracted_formula_series = maven_df[cons.FORMULA_COL].apply(get_extracted_element)
    extracted_formula_series.apply(lambda x: element_dict.update(x))

    return element_dict.keys()


def melt_df(df1):
    """
    This function melts the dataframe in long form based on some fixed columns
    Args:
        df : dataframe to be converted in long format
    Returns:
        long_form : dataframe in long format

    """
    fixed_cols = [maven_constants.NAME, maven_constants.LABEL, maven_constants.FORMULA]
    col_headers = df1.columns.tolist()
    check_column_headers(col_headers, fixed_cols)
    melt_cols = [x for x in col_headers if x not in fixed_cols]

    #try:
    
    long_form = pd.melt(df1, id_vars=fixed_cols, value_vars=melt_cols)
    
    #except KeyError():
        #raise KeyError('columns {} not found in input data'.format(','.join(fixed_cols)))
    return long_form
def maven_merge_dfs(df1, df2):
    """
    This function combines the MAVEN input file dataframe and the metadata
    file dataframe

    Args:
        input_data : MAVEN input data in form of pandas dataframe
        metadata : metadata in the form of pandas dataframe

    Returns:
        combined_data : dataframe with input data and metadata combined
    """
    long_form = melt_df(df1)
    try:
        merged_df = merge_two_dfs(long_form, df2, how='left',
                            left_on=cons.VAR_COL, right_on=maven_constants.SAMPLE)
    except KeyError:
        raise KeyError(maven_constants.SAMPLE + ' column not found in metadata')

    df_std_form = column_manipulation(merged_df)
    return df_std_form

def get_merge_df(maven_df, metadata_df):
    """
    This function merge the metadata_df with maven_df. If metadata_df
    is not present it converts the maven df in long format (standard format)
    :param mave_df: df of raw maven input file
    :param metadata_df: df of metadata info file
    :return merge_df: std wide form df
    """
    if check_df_empty(metadata_df):
        return convert_inputdata_to_stdfrom(maven_df)
    else:
        return maven_merge_dfs(maven_df, metadata_df)



def read_maven_file(maven_file_path, maven_df, metadata_df,validation_logs,summary):
    """
    This function reads maven and metadata file, convert it to df and
    checks for validation of files. If validation does not raise any
    error it returns mergedf with logs and iso-tracer data.
    :param maven_file_path: absolute path of maven raw file
    :param maven_sample_metadata_path: absolute path of metadatafile
    :return: mergedf : merge df of Maven and Metadata File
             logs: dictionary of errors and warnings
             iso-tracer : dictionary of iso-tracer details
    """
    if not check_error_present(validation_logs):
        isotracer_dict = get_isotracer_dict(maven_df)
        merged_df = get_merge_df(maven_df, metadata_df)
        unique_element_list = get_element_list(maven_df)
        raw_df_with_no_header = read_file(maven_file_path, head=None)
        list_of_columns = list(raw_df_with_no_header.iloc[0])
        duplicate_column_list = check_duplicates_in_list(list_of_columns)
        if duplicate_column_list: 
            duplicate_message = "Found column " + ','.join(duplicate_column_list) + " are duplicated"
            action_mesaage = "Column are renamed and appended"
            validation_logs['warnings']['message'].append(duplicate_message)
            validation_logs['warnings']['action'].append(action_mesaage)
        return merged_df, validation_logs, isotracer_dict, unique_element_list, summary
    else:
        return corrected_maven_df, validation_logs, None, None, summary


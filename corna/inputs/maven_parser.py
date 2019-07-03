import pandas as pd
from collections import namedtuple

from corna import helpers as hlp
from column_conventions import maven as maven_constants
import corna.constants as cons

MavenKey = namedtuple('MavenKey', 'name formula')

def get_isotope_columns_frm_label_col(df, iso_tracers):
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
            for i in hlp.chemformula_schema.parseString(formula))
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


def filter_required_col_and_get_formula_dict(df, metabolite, isotracers, required_column):
    """
    This function returns formula dict for the formula of metabolite.
    It also filters the dataframe according to columns required.
    formula- C5H4N3
    formula_dict- {'C':5,'H':4, 'N':3}

    Parameters:
        df: Dataframe to be processed
        metabolit: metabolite name
        isotracer: list of isotracers
        required_column: columns required for processing.

    Returns:
        required_df: df required for processing.
        formula: formula of metabolite
        formula_dict: dict of form- element:no. of atoms present in formula
    """
    metabolite_df=df[df[cons.NAME_COL]==metabolite]
    formula = metabolite_df.Formula.unique()
    if len(formula) == 1:
        formula_dict = hlp.parse_formula(formula[0])
    else:
        lst = metabolite_df.Formula.tolist()
        formula_unique = max(lst,key=lst.count)
        metabolite_df.Formula = formula_unique
        formula_dict = hlp.parse_formula(formula_unique)
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
    duplicate_list = [x for x in given_list if given_list.count(x) > 1 ] 
    return list(set(duplicate_list))


def create_label_column_frm_isotope_columns(df, isotracers):
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


def add_name_formula_label_col(info_df, metab, formula, iso_tracers, eleme_corr):
    """
    Adds required columns back to the na corrected Dataframe 
    Required columns include : Label, Name, Formula, indistinguishable_isotope_dict

    Args:
        info_df: df to be processed
        metab: metabolite name
        formula: formula of metabolite
        iso_tracers: list of isotracers
        eleme_corr: if user selects autodetect=False, they can give a standard
                    dict of indistinguishable elements for correction.
                    eg - {'C13':['H','O']}
    Returns:
        info_df: processed dataframe.
    """
    info_df=create_label_column_frm_isotope_columns(info_df, iso_tracers)
    info_df[cons.NAME_COL]= metab
    info_df[cons.FORMULA_COL] = formula
    info_df[cons.INDIS_ISOTOPE_COL] = str(eleme_corr)
    return info_df


def save_original_label_and_processed_label(df, iso_tracers):
    """
    This function modifies Label column from C13-Label-1 to C13N15-Label-1-0, save it as well.
    And then returns back the dataframe.

    Args:
        df: dataframe to be processed.
        iso_tracers: list of isotracers present.

    Returns:
        df: processed dataframe
    """
    df= df.filter([cons.NAME_COL, cons.FORMULA_COL, cons.LABEL_COL, cons.SAMPLE_COL, cons.INTENSITY_COL])
    df[cons.ORIGINAL_LABEL_COL]=df[cons.LABEL_COL]
    df= get_isotope_columns_frm_label_col(df, iso_tracers)
    df= create_label_column_frm_isotope_columns(df, iso_tracers)
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
    return hlp.parse_formula(formula)


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
    hlp.check_column_headers(col_headers, fixed_cols)
    melt_cols = [x for x in col_headers if x not in fixed_cols]

    try:    
        long_form = pd.melt(df1, id_vars=fixed_cols, value_vars=melt_cols)
    
    except KeyError():
        raise KeyError('columns {} not found in input data'.format(','.join(fixed_cols)))
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
        merged_df = hlp.merge_two_dfs(long_form, df2, how='left',
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



def read_maven_file(maven_file_path, maven_df, metadata_df):
    """
    This function reads maven and metadata df. If validation does not 
    raise any error it returns mergedf with logs and iso-tracer data.
    :param maven_file_path: absolute path of maven raw file
    :param maven_df: maven raw dataframe
    :param metadata_df: metadata dataframe if present otherwise empty df.
    :param validation_logs: logs after validation of input files.
    :return: mergedf : merge df of Maven and Metadata File
             isotracer_dict : dictionary of iso-tracer details
             unique_element_lis: list of unique elements present 
                                in formula column
    """
    isotracer_dict = get_isotracer_dict(maven_df)
    merged_df = get_merge_df(maven_df, metadata_df)
    unique_element_list = get_element_list(maven_df)
    return merged_df, isotracer_dict, unique_element_list


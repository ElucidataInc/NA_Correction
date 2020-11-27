"""
Note:
    NA values are referenced from:
    CRC Handbook of Chemistry and Physics (83rd ed.). Boca Raton, FL: CRC Press. ISBN 0-8493-0483-0.
    by Lide, D. R.
    ed. (2002)
"""

import os
import json

elementdata = os.path.join(os.path.dirname(__file__),'element_data.json')
with open(elementdata) as data_file:
    data = json.load(data_file)

ELE_ATOMIC_WEIGHTS = data["element_mol_weight_dict"]
ISOTOPE_NA_MASS = data["isotope_na_mass"]
UPPER_CASE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
LOWER_CASE = "abcdefghijklmnopqrstuvwxyz"
DIGITS = "0123456789"
VAR_COL = "variable"
VAL_COL = "value"
ELEMENT_SYMBOL = ('[A-Z][a-z]?')
LEVEL_0_COL = "level_0"
LEVEL_1_COL = "level_1"
INTENSITY_COL = 'Intensity'
SAMPLE_COL = 'Sample'
ORIGINAL_LABEL_COL= 'Original_label'
NA_CORRECTED_COL= 'NA Corrected'
NA_CORRECTED_WITH_ZERO= 'NA Corrected with zero'
FRACTIONAL_ENRICH= 'Fractional enrichment'
SAMPLE_NAME = 'Sample Name'
PARENT_COL = 'Unlabeled Fragment'
LABEL_COL = 'Label'
FRAG_COL = 'frag_keys'
MASSINFO_COL = 'Mass Info'
ISOTRACER_COL = 'Isotopic Tracer'
FORMULA_COL = 'Formula'
VALID_STATE = 'correct'
MISSING_STATE = 'missing'
DUPLICATE_STATE = 'duplicate'
LABEL_STATE_INVALID = 'invalid label format'
LABEL_STATE_NOT_CORRECT = 'a required format "C13-Label-1"'
LABEL_STATE_NOT_FORMULA = 'no label in formula'
LABEL_STATE_NUMBER_MORE_FORMULA = 'element_in_label_more_than_formula'
ELEMENT_LIST = set(['C','N','H','S', 'O'])
UNLABELLED_LABEL = 'C12 PARENT'
FORMULA_STATE_INVALID = 'incorrect formula'
COLUMN_STATE = 'state'
COLUMN_ROW = 'row_number'
COLUMN_NAME = 'column_name'
COMPONENT_NAME = 'Component Name'
NAME_COL = 'Name'
PARENT_MASS_MOL= 'Parent_mass'
DAUGHTER_MASS_MOL= 'Daughter_mass'
PARENT_NUM_ATOMS= 'Number of atoms parent'
DAUGHTER_NUM_ATOMS= 'Number of atoms daughter'
PARENT_NUM_LABELED_ATOMS= 'Number of labeled atoms parent'
DAUGHTER_NUM_LABELED_ATOMS= 'Number of labeled atoms daughter'
PARENT_MASS_ISO= 'isotopic_mass_parent'
DAUGHTER_MASS_ISO= 'isotopic_mass_daughter'
UNLABELED= 'Unlabeled'
BACKGROUND_CORRECTED= 'Background Corrected'
BACKGROUND_WITH_ZERO= 'Background Corrected with zero'
PARENT_FORMULA_COL= 'Parent_formula'
NA_LCMS = 'na_lcms'
NA_MSMS = 'na_msms'
UNLABELLED_LABEL_DICT = {'C': 0,'N': 0}
INTENSITY_STATE_NEGATIVE = 'negative intensity'
INTENSITY_STATE_INVALID = 'invalid_intensity_value'
VALIDATION_WARNING = 'warnings'
VALIDATION_ERROR = 'errors'
VALIDATION_MESSAGE = 'message'
VALIDATION_ACTION = 'action'
VALIDATION_ACTION_DROP = 'DROP'
VALIDATION_ACTION_FILL_NA = 'FILL_NA'
VALIDATION_ACTION_STOP = 'Stop_Tool'
VALIDATION_ACTION_OK = 'All_Ok'
VALIDATION_ACTION_ROW_WISE = 'Row_Wise_Action'
WARNING_STATE = [DUPLICATE_STATE,MISSING_STATE]
VALIDATION_COLUMN_NAME = 'column'
VALIDATION_MSG_ROW_DROPPED = "Row is Dropped"
VALIDATION_MSG_FILL_NA = "Missing value of columns replaced with 0"
MOL_MASS_VALIDATE = 'Molecular weight of a metabolite cannot be zero'
PPM_REQUIREMENT_VALIDATION = 'The ppm requirement is at the boderline for '
MISSING_COMPONENTS = "missing_components"
##Keys in ISOTOPE_NA_MASS dictionary
KEY_NA = "naValue"
KEY_AMU = "amu"
KEY_ELE = "element"
KEY_NAT_ISO = "naturalIsotope"
KEY_AMU_DIFF = 'amu_diff'
ORIGINAL_FILENAME = "Original Filename"
BACKGROUND_SAMPLE = "Background Sample"
FORMULA_COL_METADATAFILE = ['Formula', 'Parent Formula']
AREA_COLUMN_RAWFILE = ['Area']
COLUMN_ISOTOPE_TRACER = 'Isotopic Tracer'
ISOTOPE_VALUES = ['C13', 'N15', 'H2', 'S34']
BORDERLINE_LIMIT = 0.5
PATTERN_MASSINFO_COL = '\d+.*\d+\s*\/\s*\d+.*\d+'
SAMPLE_NAME_COL_PATTERN = '^([a-zA-Z0-9_\s\-]*)$'
## Dict storing mass diff between isotopes
#formula used = {(O17-O16) - (C13-C12)}
MASS_DIFF_DICT= {
  "O17": {
      "C": 0.0008622426,
      "N": 0.007182187,
      "H": 0.0020596653,
      "S34": 0.0126382609,
      "Si29": 0.0046489129,
      "Si30": 0.0115905234,
      "S33": 0.0048293204
    },
      "H2": {
      "O17": 0.0020596653,
      "C": 0.0029219079,
      "N": 0.0092418523,
      "S34": 0.0167575915,
      "Si29": 0.0067085782,
      "Si30": 0.015709854,
      "S33": 0.0068889857,
      "O18": 0.008307111
    },
         "N15": {
      "O17": 0.007182187,
      "C": 0.0063199444,
      "H": 0.0092418523,
      "S34": 0.0017261132,
      "Si29": 0.0025332741,
      "Si30": 0.0027738507,
      "S33": 0.0023528666,
      "O18": 0.0101765936
    },
     "C13": {
      "O17": 0.0008622426,
      "N": 0.0063199444,
      "H": 0.0029219079,
      "S34": 0.0109137756,
      "Si29": 0.0037866703,
      "Si30": 0.0098660381,
      "S33": 0.0039670778,
      "O18": 0.0024632952
    },
    "S34": {
      "O17": 0.0126382609,
      "C": 0.0109137756,
      "N": 0.0017261132,
      "H": 0.0167575915,
      "Si29": 0.003340435,
      "Si30": 0.0010477375,
      "O18": 0.0084504804
    },
     "Si29": {
      "O17": 0.0046489129,
      "C": 0.0037866703,
      "N": 0.0025332741,
      "H": 0.0067085782,
      "S34": 0.003340435,
      "S33": 0.0001804075,
      "O18": 0.0051100454
    },
     "Si30": {
      "O17": 0.0115905234,
      "C": 0.0098660381,
      "N": 0.0027738507,
      "H": 0.015709854,
      "S34": 0.0010477375,
      "S33": 0.0019318825,
      "O18": 0.0074027429
    },
     "S33": {
      "O17": 0.0048293204,
      "C": 0.0039670778,
      "N": 0.0023528666,
      "H": 0.0068889857,
      "Si29": 0.0001804075,
      "Si30": 0.0019318825,
      "O18": 0.0054708604
    },
     "O18": {
      "C": 0.0024632952,
      "N": 0.0101765936,
      "H": 0.008307111,
      "S34": 0.0084504804,
      "Si29": 0.0051100454,
      "Si30": 0.0074027429,
      "S33": 0.0054708604
    }
}
## Isotope dictionary
ISOTOPE_DICT = {'O': ['O17', 'O18'], 'S': ['S33', 'S34'], 'Si': ['Si29', 'Si30']}

##List for isotopic elements that are added in na_dict
ISOTOPE_LIST = ['O17', 'O18', 'S33', 'S34', 'Si29', 'Si30']

##Output column names
INDIS_ISOTOPE_COL = 'Indistinguishable_isotope'
POOL_TOTAL_COL = 'Pool_total'
METABOLITE_NAME = 'metab_name'

##summary tab
SUMMARY_LABEL = 'label'
SUMMARY_VAL = 'value'
SUMMARY_TITLE = 'title'
SUMMARY = 'summary'
RAW_FIELD_SUMMARY_LIST = ['Number of rows', 'Number of samples', 'Number of cohorts', 'Number of metabolites']
META_FIELD_SUMMARY_LIST = ['Number of fragments', 'Number of unlabeled fragments', 'isotopic tracer']
SAMPLE_FIELD_SUMMARY_LIST = ['Number of background samples', 'Fields in metadata']
LCMS_RAW_FIELD_SUMMARY = ['Number of metabolites', 'Number of samples', 'Number of blank intensity cells', 'Number of rows']
LCMS_META_FILED_SUMMARY = ['Fields in metadata', 'Number of rows in metadata']

## file type
RAW_MSMS = 'InputData'
META_MSMS = 'MetaData'
SMP_MSMS = 'SampleData'
RAW_LCMS = 'Input_Data'
META_LCMS = 'Meta_Data'

FILE_PATH = 'file_path'
FORMULA = 'Formula'
FORMULA_COL_PATTERN = 'C+'

SAMPLE_METADATA_REQUIRED_COLS = ['Original Filename', 'Sample Name', 'Background Sample']

METADATA_MQ_REQUIRED_COLS = ['Component Name', 'Unlabeled Fragment', 'Isotopic Tracer', 'Formula',
                             'Parent Formula']

RAW_FILE_REQUIRED_COLS_WITH_SAMPLE_METADATA = ['Original Filename', 'Area',
                             'Mass Info', 'Sample Name', 'Component Name']
RAW_FILE_REQUIRED_COLS_WITHOUT_SAMPLE_METADATA = ['Original Filename', 
                             'Area', 'Mass Info', 'Component Name']

RAW_MQ_DICT = {
    'file_path': None,
    'required_columns': [],
    'warnings': {
        'missing': 'FILL_NA',
        'duplicate': 'DROP',
    },
    'functions': {
        'numerical': {'column_list': AREA_COLUMN_RAWFILE,
                      'negative state': 'negative intensity',
                      'invalid state': 'invalid num'},
        'pattern_match': {'column_name': MASSINFO_COL,
                          'regex_pattern': PATTERN_MASSINFO_COL,
                          'state':'not in correct format'},
        'missing_data': {'state': 'missing'},
    }
}

METADATA_MQ_DICT = {
    'file_path': '',
    'required_columns': METADATA_MQ_REQUIRED_COLS,
    'warnings': {
        'missing': 'FILL_NA',
        'duplicate': 'DROP',
        'pattern_match': {'column_name': FORMULA,
                          'regex_pattern': FORMULA_COL_PATTERN,
                          'state': 'not in correct format'},
    },
    'functions': {
        'chemical_formula': {'column_list': FORMULA_COL_METADATAFILE,
                            'state': 'invalid formula'},
        'value_in_constant': {'column_name': ISOTRACER_COL,
                              'constant_list' : ISOTOPE_VALUES,
                              'state' : 'invalid'},
        'missing_data': {'state': 'missing'},
        }
}

SAMPLE_METADATA_DICT = {
    'file_path': None,
    'required_columns': SAMPLE_METADATA_REQUIRED_COLS,
    'warnings': {
        'missing': 'FILL_NA',
        'duplicate': 'DROP',
    },
    'functions': {
        'missing_data': {'state': 'missing'},
        }
}

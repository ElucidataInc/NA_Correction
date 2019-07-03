# TODO : many functions need to be documneted refer issue NCT-303

import os
import warnings

import pandas as pd

from corna.inputs.column_conventions import multiquant
from corna import constants as const
import corna.helpers as hlp

from collections import namedtuple
import os
import warnings

from datum import algorithms as dat_alg
import pandas as pd

from .column_conventions import multiquant
from corna import constants
from corna import summary as sm
from ..constants import INTENSITY_COL
from corna.inputs import validation
from ..data_model import standard_model
from ..helpers import read_file, get_unique_values, check_column_headers
from ..isotopomer import bulk_insert_data_to_fragment

Multiquantkey = namedtuple('MultiquantKey', 'name formula parent parent_formula')
validated_raw_tuple = namedtuple('validated_raw_mq', 'df logs')
validated_metadata_tuple = namedtuple('validated_metadata_mq', 'df logs')
metadata_mq_tuple = namedtuple('metadata_mq', 'df logs')

def add_mass_and_no_of_atoms_info_frm_label(df):
    """
    This function adds column of parent isotopic mass, daughter isotopic mass, parent molecular mass,
    daughter molecular mass, no of labeled atoms parent and daughter,
    total no of atoms parent and daughter and isotracer column to the dataframe.
    """
    msms_df = get_parent_daughter_iso_mass_col_and_isotracer_from_label(df)
    isotracer= msms_df[multiquant.ISOTRACER].unique()
    iso_elem= hlp.get_isotope_element(isotracer[0])
    input_df=pd.DataFrame()
    final_df=pd.DataFrame()

    input_df= get_mol_mass_and_number_of_atoms(msms_df, multiquant.PARENT_FORM, 'PARENT', isotracer[0],iso_elem, input_df)
    final_df= get_mol_mass_and_number_of_atoms(input_df, multiquant.FORMULA, 'DAUGHTER', isotracer[0],iso_elem, final_df)

    return final_df, isotracer


def get_num_labeled_atoms(isotope, isotopic_mass, molecular_mass):
    """
    This function returns no of labeled atoms of isotracer from isotopic mass and molecular mass of compound.

    Parameters:
        isotope: isotope element used for labeling
        isotopic_mass: mass of compound due to presence of isotope.
        molecular_mass: molecular mass of compound.

    Returns:
        number_label: number of labeled atoms in the formula.

    """
    nat_iso = hlp.get_isotope_natural(isotope)
    if nat_iso == isotope:
        return 0
    atom_excess_mass = hlp.get_isotope_mass(
                isotope) - hlp.get_isotope_mass(nat_iso)
    label= (isotopic_mass - molecular_mass) / atom_excess_mass
    label.round(2)
    number_label= label.astype(int)
    return number_label


def get_parent_daughter_iso_mass_col_and_isotracer_from_label(msms_df):
    """
    This function creates columns of parent isotopic mass, daughter isotopic mass and isotracer from 
    the label column.
    """
    s=msms_df[multiquant.LABEL].apply(lambda x: x.split('_'))
    msms_df[const.PARENT_MASS_ISO] = pd.to_numeric(s.apply(lambda x: x[1]), errors='coerce')
    msms_df[const.DAUGHTER_MASS_ISO] = pd.to_numeric(s.apply(lambda x: x[2]), errors='coerce')
    msms_df[multiquant.ISOTRACER]= s.apply(lambda x: x[0])   
    msms_df.rename(columns={multiquant.PARENT_FORMULA: multiquant.PARENT_FORM}, inplace=True)
    return msms_df


def get_mol_mass_and_number_of_atoms(df, formula_col, formula_info, iso_tracer,iso_elem, out_df):
    """
    This function creates the columns of molecular mass of parent and daughter, total no of atoms of
    isotope in parent formula and daughter formula, number of labeled atoms of isotope in parent 
    and daughter formula.

    Args:
        df: input df
        formula_col: column name of formula column
        formula_info: PARENT formula or DAUGHTER formula
        iso_tracer: isotracer present in formula
        iso_elem: element name of isotracer

    Returns:
        out_df: output dataframe
    """
    if formula_info=='PARENT':
        mol_mass_col= const.PARENT_MASS_MOL
        iso_mass_col= const.PARENT_MASS_ISO
        num_atom_col= const.PARENT_NUM_ATOMS
        num_labeled_atoms_col= const.PARENT_NUM_LABELED_ATOMS
    else:
        mol_mass_col= const.DAUGHTER_MASS_MOL
        iso_mass_col= const.DAUGHTER_MASS_ISO
        num_atom_col= const.DAUGHTER_NUM_ATOMS
        num_labeled_atoms_col= const.DAUGHTER_NUM_LABELED_ATOMS

    for formula in df[formula_col].unique():
        form_df= df[df[formula_col]==formula]
        form_dict= hlp.parse_formula(formula)
        try:
            form_df.loc[:,num_atom_col] = form_dict[iso_elem]
        except KeyError:
            form_df.loc[:,num_atom_col] = 0
        mass= hlp.get_mol_weight(formula)   
        form_df.loc[:,mol_mass_col] = mass
        form_df.loc[:,mol_mass_col]= form_df[mol_mass_col].round().astype(int)
        out_df= out_df.append(form_df)

    out_df[num_labeled_atoms_col]= get_num_labeled_atoms(iso_tracer, out_df[iso_mass_col], out_df[mol_mass_col])
    return out_df


def mq_merge_meta(input_data, metadata):
    """
    This function combines the MQ input file dataframe and the metadata
    file dataframe

    Args:
        input_data : MQ input data in form of pandas dataframe

        metadata : metadata in the form of pandas dataframe

    Returns:
        merged_df : dataframe with input data and metadata combined
    """

    try:
        merged_df = input_data.merge(metadata, how='inner',
                                     left_on=multiquant.MQ_FRAGMENT,
                                     right_on=multiquant.MQ_FRAGMENT)
        if merged_df.empty:
            raise Exception('Empty Merge, no common entries to process, please check input files')
    except KeyError:
        raise KeyError('Missing columns: ' + multiquant.MQ_FRAGMENT)
    merged_df[multiquant.MASSINFO] = merged_df[multiquant.MASSINFO].str.replace(' / ', "_")
    return merged_df


def merge_samples(merged_df, sample_metadata):
    """
    This function merge the raw input dataframe with the sample metadata 
    dataframe.
    """
    if sample_metadata is not None:
        col_headers_sample = sample_metadata.columns.values
        col_headers_merged = merged_df.columns.values
        bg_corr_col_names_sample = [multiquant.BACKGROUND, multiquant.MQ_COHORT_NAME]
        bg_corr_col_names_merged = [multiquant.MQ_COHORT_NAME]
        try:
            hlp.check_column_headers(col_headers_sample, bg_corr_col_names_sample)
            hlp.check_column_headers(col_headers_merged, bg_corr_col_names_merged)
            assert set(sample_metadata[multiquant.BACKGROUND]).issubset(set(sample_metadata[multiquant.MQ_SAMPLE_NAME]))
            merged_df = merged_df.merge(sample_metadata, how='inner',
                                        on=[multiquant.MQ_SAMPLE_NAME, multiquant.MQ_COHORT_NAME])
            if merged_df.empty:
                raise Exception('Empty Merge, no common entries to process, please check input files')
        except AssertionError:
            warnings.warn("Background Correction can't be performed")
            merged_df = merged_df.merge(sample_metadata, how='inner',
                                        on=[multiquant.MQ_SAMPLE_NAME])

    # first change Sample Name to Cohort Name, then Original Filename to Sample
    # refer to multiquant raw output
    merged_df.rename(
        columns={multiquant.MQ_COHORT_NAME: multiquant.COHORT}, inplace=True)
    merged_df.rename(
        columns={multiquant.MQ_SAMPLE_NAME: multiquant.SAMPLE}, inplace=True)

    remove_stds = remove_mq_stds(merged_df)
    remove_stds.rename(columns={"Area": multiquant.INTENSITY}, inplace=True)
    return remove_stds

def get_background_samples(sample_metadata, sample_name, background_sample):
    sample_background = sample_metadata.set_index(sample_name).to_dict()
    return sample_background[background_sample]

def mq_merge_dfs(input_data, metadata, sample_metadata):
    """
    First merge metadata with the raw input file and then merge sample metadata file.
    """
    merged_data = mq_merge_meta(input_data, metadata)
    return merge_samples(merged_data, sample_metadata)


def remove_mq_stds(merged_df):
    """
    This function removes the standard samples from multiquant data
    """
    try:
        merged_df = merged_df[not merged_df[multiquant.COHORT].str.contains("std")]
    except:
        warnings.warn('Std samples not found in' + multiquant.COHORT + ' column')

    merged_df[multiquant.LABEL] = merged_df[multiquant.ISOTRACER] + "_" + merged_df[multiquant.MASSINFO]
    merged_df.rename(
        columns={multiquant.PARENT: multiquant.NAME}, inplace=True)
    merged_df.pop(multiquant.MASSINFO)
    merged_df.pop(multiquant.ISOTRACER)
    return merged_df


def get_replicates(sample_metadata, sample_name, cohort_name, background_sample):
    """
    This function returns the list of sample names which belong
    to the same cohort.

    Parameters
    ----------
        sample_metadata : DataFrame
            sample metadata df
        sample_name : string
            column name of column containing information of original filename
        cohort name : string
            column name of column containing information of cohort name.
        background sample: string
            column name of column containing information of corresponding background file names.
    
    Returns
    -------
        replicate_groups: list of sample names which belong to the same cohort.
    """
    sample_index_df = sample_metadata.set_index(sample_name)
    sample_index_df['Background Cohort'] = sample_index_df[
        background_sample].map(sample_index_df[cohort_name])
    replicate_groups = []
    # std should not be present in sample_metadata
    cohort_list = hlp.get_unique_values(sample_index_df, 'Background Cohort')
    for cohorts in cohort_list:
        newdf = sample_index_df[sample_index_df[
                                    'Background Cohort'] == cohorts]
        replicate_groups.append(hlp.get_unique_values(newdf, background_sample))
    return replicate_groups

def frag_key(df):
    """
    This function creates a fragment key column in merged data based on parent information.
    """

    def _extract_keys(x):
        return Multiquantkey(x[multiquant.MQ_FRAGMENT],
                             x[multiquant.FORMULA],
                             x[multiquant.NAME],
                             x[multiquant.PARENT_FORMULA])

    try:
        df[multiquant.FRAG] = df.apply(_extract_keys, axis=1)
    except KeyError:
        raise KeyError('Missing columns in data')
    return df

def merge_mq_metadata(mq_df, metdata, sample_metdata):
    """
    This function merge the raw input dataframe with the metadata and 
    sample metadata.
    Args:
        mq_df: raw input df
        metdata: metadata df
        sample_metdata: sample metadata df

    Returns:
        merged_df: merged dataframe
        list_of_replicates: list of sample names belonging to same cohort.
    """
    merged_data = mq_merge_dfs(mq_df, metdata, sample_metdata)
    merged_data.fillna(0, inplace=True)
    list_of_replicates = []
    sample_background = []

    if sample_metdata is not None:
        col_headers = merged_data.columns.values
        bg_corr_col_names = [multiquant.BACKGROUND, multiquant.COHORT]
        try:
            hlp.check_column_headers(col_headers, bg_corr_col_names)
        except AssertionError:
            return merged_data, list_of_replicates, sample_background
        # consider only those sample names which are present in raw file
        sample_metdata = sample_metdata[sample_metdata[multiquant.MQ_SAMPLE_NAME].isin(merged_data[multiquant.SAMPLE])]
        list_of_replicates = get_replicates(
            sample_metdata, multiquant.MQ_SAMPLE_NAME, multiquant.MQ_COHORT_NAME, multiquant.BACKGROUND)
        sample_background = get_background_samples(
            sample_metdata, multiquant.MQ_SAMPLE_NAME, multiquant.BACKGROUND)
    return merged_data, list_of_replicates, sample_background


def mq_df_to_fragmentdict(merged_df, intensity_col=INTENSITY_COL):
    frag_key_df = frag_key(merged_df)
    print("frag_key dataframe")
    #print(frag_key_df.head())
    std_model_mq = standard_model(frag_key_df, intensity_col)
    #print(std_model_mq)
    metabolite_frag_dict = {}
    for frag_name, label_dict in std_model_mq.iteritems():
        #print(frag_name, label_dict)
        curr_frag_name = Multiquantkey(frag_name.name, frag_name.formula,
                                       frag_name.parent, frag_name.parent_formula)
        #print(curr_frag_name)
        if metabolite_frag_dict.has_key(frag_name.parent):
            print("it has key")
            metabolite_frag_dict[frag_name.parent].update(bulk_insert_data_to_fragment(curr_frag_name,
                                                                                       label_dict, mass=True))
            #print("metabolite frag dict printed")
            #rint(metabolite_frag_dict)
        else:
            print("no key")
            metabolite_frag_dict[frag_name.parent] = bulk_insert_data_to_fragment(curr_frag_name,
                                                                                  label_dict, mass=True)
    return metabolite_frag_dict

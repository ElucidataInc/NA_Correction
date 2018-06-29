from collections import namedtuple
import pandas as pd
import numpy as np
import re

import corna.constants as const
import corna.helpers as hlp
from corna.inputs.column_conventions import multiquant 

def get_num_labeled_atoms(isotope, isotopic_mass, molecular_mass):
    """
    This function returns no of labeled atoms of isotracer from isotopic mass and molecular mass of compound.
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

def get_input_format(msms_df):
    """
    This function creates columns of parent isotopic mass, daughter isotopic mass, isotracer from 
    the label column.
    """
    s=msms_df[multiquant.LABEL].apply(lambda x: x.split('_'))
    msms_df[const.PARENT_MASS_ISO] = pd.to_numeric(s.apply(lambda x: x[1]), errors='coerce')
    msms_df[const.DAUGHTER_MASS_ISO] = pd.to_numeric(s.apply(lambda x: x[2]), errors='coerce')
    msms_df[multiquant.ISOTRACER]= s.apply(lambda x: x[0])   
    msms_df[const.NA_CORRECTED_COL]=0.0
    msms_df.rename(columns={multiquant.PARENT_FORMULA: multiquant.PARENT_FORM}, inplace=True)
    return msms_df

def get_mass_and_number_of_atoms(df, formula_col, formula_info, iso_tracer,iso_elem, out_df):
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


def na_correction_mimosa(msms_df):
    """
    This function performs na correction on the input data frame for LCMS/MS file.
    Args:
        msms_df: Dataframe which contains intensities to be corrected.

    Returns:
        output_df: na corrected dataframe
        
    """
    REQUIRED_COL= [multiquant.FORMULA, multiquant.LABEL, multiquant.NAME, multiquant.SAMPLE, multiquant.COHORT,
                     multiquant.MQ_FRAGMENT, multiquant.INTENSITY, multiquant.PARENT_FORM,const.NA_CORRECTED_COL]
    msms_df= get_input_format(msms_df)
    isotracer= msms_df[multiquant.ISOTRACER].unique()
    iso_elem= hlp.get_isotope_element(isotracer[0])
    input_df=pd.DataFrame()
    final_df=pd.DataFrame()
    output_df= pd.DataFrame()
    metab_dict={}

    na= hlp.get_isotope_na(isotracer[0])

    input_df= get_mass_and_number_of_atoms(msms_df, multiquant.PARENT_FORM, 'PARENT', isotracer[0],iso_elem, input_df)
    final_df= get_mass_and_number_of_atoms(input_df, multiquant.FORMULA, 'DAUGHTER', isotracer[0],iso_elem, final_df)

    final_df['A']=(1 + na * (final_df[const.PARENT_NUM_ATOMS]-final_df[const.PARENT_NUM_LABELED_ATOMS]))
    final_df['B']= na * ((final_df[const.PARENT_NUM_ATOMS]-final_df[const.DAUGHTER_NUM_ATOMS]) -\
                         (final_df[const.PARENT_NUM_LABELED_ATOMS]-final_df[const.DAUGHTER_NUM_LABELED_ATOMS]-1))
    final_df['C']=  na * (final_df[const.DAUGHTER_NUM_ATOMS]-final_df[const.DAUGHTER_NUM_LABELED_ATOMS]+1)

    final_df.drop([const.PARENT_MASS_MOL, const.DAUGHTER_MASS_MOL, const.PARENT_NUM_ATOMS,
                const.DAUGHTER_NUM_ATOMS, const.DAUGHTER_NUM_LABELED_ATOMS, const.PARENT_NUM_LABELED_ATOMS], axis=1, inplace=True)

    for samp in final_df.Sample.unique():
        """
        Create metabolite dictionary of the form:
            {'SAMPLE 2_10':{
                (191, 111): 2345.75, (192, 111):5644.847
                }
            }
        """
        metab_df = final_df[final_df[multiquant.SAMPLE]==samp]
        frag_dict={}
        for index, row in metab_df.iterrows():
            frag_dict[(row[const.PARENT_MASS_ISO],row[const.DAUGHTER_MASS_ISO])]=row[multiquant.INTENSITY]
    
        metab_dict[samp]= frag_dict

    for samp in final_df.Sample.unique():
        metab_df = final_df[final_df[multiquant.SAMPLE]==samp]
        frag= metab_dict[samp]  
     
        for index, row in metab_df.iterrows():
            m_n= row[const.DAUGHTER_MASS_ISO]
            m_1_n= row[const.PARENT_MASS_ISO]-1
            m_n_1= row[const.DAUGHTER_MASS_ISO]-1
            intensity_m_n= row[const.INTENSITY_COL]
            try:
                intensity_m_1_n= frag[m_1_n, m_n]
            except KeyError:
                intensity_m_1_n=0
            try:
                intensity_m_1_n_1= frag[m_1_n, m_n_1]
            except KeyError:
                intensity_m_1_n_1= 0
        
            corrected= intensity_m_n * row['A']  - intensity_m_1_n * row['B'] -\
                                                intensity_m_1_n_1 * row['C']
            metab_df.set_value(index=index, col=const.NA_CORRECTED_COL, value= corrected)

        output_df=output_df.append(metab_df) 
    output_df= output_df.filter(REQUIRED_COL)

    return output_df



  


        






    
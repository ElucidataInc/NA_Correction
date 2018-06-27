from collections import namedtuple
import pandas as pd
import numpy as np
import re

import corna.constants as const
import corna.helpers as hlp

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
    s=msms_df['Label'].apply(lambda x: x.split('_'))
    msms_df['isotopic_mass_parent'] = pd.to_numeric(s.apply(lambda x: x[1]), errors='coerce')
    msms_df['isotopic_mass_daughter'] = pd.to_numeric(s.apply(lambda x: x[2]), errors='coerce')
    msms_df['isotracer']= s.apply(lambda x: x[0])   
    msms_df['NA Corrected']=0.0
    msms_df.rename(columns={'Parent Formula': 'Parent_formula'}, inplace=True)
    return msms_df


def na_correction_mimosa(msms_df):
    """
    This function performs na correction on the input data frame for LCMS/MS file.
    Args:
        msms_df: Dataframe which contains intensities to be corrected.

    Returns:
        output_df: na corrected dataframe
        
    """
    msms_df= get_input_format(msms_df)
    required_columns=['Formula','Label','Name', 'Sample', 'Cohort Name', 'Component Name', 'Intensity', 'NA Corrected', 'Parent_formula']
    isotracer= msms_df.isotracer.unique()
    input_df=pd.DataFrame()
    final_df=pd.DataFrame()
    output_df= pd.DataFrame()
    metab_dict={}

    na= hlp.get_isotope_na(isotracer[0])

    #create column 'p' and 'Parent_mass' which contains no. of atoms of isotracer and molecular mass of 
    #parent respectively. 
    for formula in msms_df.Parent_formula.unique():
        df= msms_df[msms_df['Parent_formula']==formula]
        form_dict= hlp.parse_formula(formula)
        try:
            df.loc[:,'p'] = form_dict['C']
        except KeyError:
            df.loc[:,'p'] = 0
        mass= hlp.get_mol_weight(formula)   
        df.loc[:,'Parent_mass'] = mass
        df.loc[:,'Parent_mass']= df['Parent_mass'].round().astype(int)
        input_df= input_df.append(df)

    #create column 'd' and 'Daughter_mass' which contains no. of atoms of isotracer and molecular mass of 
    #daughter respectively. 
    for formula in input_df.Formula.unique():
        df1=input_df[input_df['Formula']==formula]
        form_dict= hlp.parse_formula(formula)
        try:
            df1.loc[:,'d'] = form_dict['C']
        except KeyError:
            df1.loc[:,'d'] = 0
        mass= hlp.get_mol_weight(formula)
        df1.loc[:,'Daughter_mass'] = mass
        df1.loc[:,'Daughter_mass']= df1['Daughter_mass'].round().astype(int)
        final_df= final_df.append(df1) 

    final_df['m']= get_num_labeled_atoms('C13', final_df['isotopic_mass_parent'], final_df['Parent_mass'])
    final_df['n']= get_num_labeled_atoms('C13', final_df['isotopic_mass_daughter'], final_df['Daughter_mass'])

    final_df['A']=(1 + na * (final_df['p']-final_df['m']))
    final_df['B']= na * ((final_df['p']-final_df['d']) - (final_df['m']-final_df['n']-1))
    final_df['C']=  na * (final_df['d']-final_df['n']+1)

    final_df.drop(['Parent_mass', 'Daughter_mass', 'p','d', 'm', 'n'], axis=1, inplace=True)

    for samp in final_df.Sample.unique():
        """
        Create metabolite dictionary of the form:
            {'SAMPLE 2_10':{
                (191, 111): 2345.75, (192, 111):5644.847
                }
            }
        """
        metab_df = final_df[final_df['Sample']==samp]
        frag_dict={}
        for index, row in metab_df.iterrows():
            frag_dict[(row['isotopic_mass_parent'],row['isotopic_mass_daughter'])]=row['Intensity']
    
        metab_dict[samp]= frag_dict

    for samp in final_df.Sample.unique():
        metab_df = final_df[final_df['Sample']==samp]
        frag= metab_dict[samp]  
     
        for index, row in metab_df.iterrows():
            m_n= row['isotopic_mass_daughter']
            m_1_n= row['isotopic_mass_parent']-1
            m_n_1= row['isotopic_mass_daughter']-1
            intensity_m_n= row['Intensity']
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
            metab_df.set_value(index=index, col='NA Corrected',value= corrected)

        output_df=output_df.append(metab_df) 
    output_df= output_df.filter(required_columns)
    
    return output_df

  


        






    
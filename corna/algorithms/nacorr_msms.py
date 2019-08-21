import pandas as pd

import corna.constants as const
from corna.helpers import get_isotope_na
from corna.inputs.column_conventions import multiquant
from corna.inputs import multiquant_parser


def na_correction_mimosa(msms_df, isBackground,
                         isotope_dict=const.ISOTOPE_NA_MASS):
    """
    This function performs na correction on the input data frame for
    LCMS/MS file. For performing na correction, it first calculates
    total number of atoms and number of labeled atoms present in both
    parent formula and daughter formula of each compound. Then for
    each sample it corrects the intensity values of each compound one
    by one using formula which includs the number of atoms, number of
    labeled atoms and intensity of M+0 isotope.
    Args:
        msms_df: Dataframe which contains intensities to be corrected.
        isBackground: boolean- True if background correction performed on
        dataframe otherwise False
        isotope_dict: Dictionary of NA values of isotopes.
                      Ex: na_dict={'naValue':{'C13':[0.9889,0.0111],
                                              'N':[0.9964,0.0036],
                                              'O':[0.9976,0.0004,0.002],
                                              'H':[0.99985,0.00015],
                                              'S':[0.95,0.0076,0.0424],
                                            }}

    Returns:
        output_df: na corrected dataframe

    """
    if isBackground:
        final_df = msms_df
        isotracer = msms_df[multiquant.ISOTRACER].unique()
        intensity_col = const.BACKGROUND_WITH_ZERO
    else:
        final_df, isotracer = \
            multiquant_parser.add_mass_and_no_of_atoms_info_frm_label(msms_df)
        intensity_col = const.INTENSITY_COL

    final_df[const.NA_CORRECTED_COL] = 0.0
    output_df = pd.DataFrame()

    na = get_isotope_na(isotracer[0], isotope_dict)

    final_df['A'] = (1 + na * (final_df[const.PARENT_NUM_ATOMS] -
                               final_df[const.PARENT_NUM_LABELED_ATOMS]))

    final_df['B'] = na * ((final_df[const.PARENT_NUM_ATOMS] -
                           final_df[const.DAUGHTER_NUM_ATOMS]) -
                          (final_df[const.PARENT_NUM_LABELED_ATOMS] -
                           final_df[const.DAUGHTER_NUM_LABELED_ATOMS] - 1))

    final_df['C'] = na * (final_df[const.DAUGHTER_NUM_ATOMS] -
                          (final_df[const.DAUGHTER_NUM_LABELED_ATOMS] - 1))

    final_df.drop([const.PARENT_MASS_MOL, const.DAUGHTER_MASS_MOL,
                   const.PARENT_NUM_ATOMS, const.DAUGHTER_NUM_ATOMS,
                   const.DAUGHTER_NUM_LABELED_ATOMS,
                   const.PARENT_NUM_LABELED_ATOMS], axis=1, inplace=True)

    metabolites = final_df[const.NAME_COL].unique()
    for metab in metabolites:
        metabolite_df = final_df[final_df[const.NAME_COL]
                                 == metab].reset_index()
        for samp in metabolite_df.Sample.unique():
            """
            Create metabolite dictionary of the form:
                {'SAMPLE 2_10':{
                    (191, 111): 2345.75, (192, 111):5644.847
                    }
                }
            """
            sample_df = metabolite_df[metabolite_df[multiquant.SAMPLE] == samp]
            frag_dict = {}
            for index, row in sample_df.iterrows():
                frag_dict[(row[const.PARENT_MASS_ISO],
                           row[const.DAUGHTER_MASS_ISO])] = row[intensity_col]
                m_n = row[const.DAUGHTER_MASS_ISO]
                m_1_n = row[const.PARENT_MASS_ISO]-1
                m_n_1 = row[const.DAUGHTER_MASS_ISO]-1
                intensity_m_n = row[intensity_col]
                try:
                    intensity_m_1_n = frag_dict[m_1_n, m_n]
                except KeyError:
                    intensity_m_1_n = 0
                try:
                    intensity_m_1_n_1 = frag_dict[m_1_n, m_n_1]
                except KeyError:
                    intensity_m_1_n_1 = 0
                corrected = intensity_m_n * row['A'] - \
                    intensity_m_1_n * row['B'] - intensity_m_1_n_1 * row['C']
                sample_df.set_value(
                    index=index, col=const.NA_CORRECTED_COL, value=corrected)
            output_df = pd.concat([output_df, sample_df])
    return output_df

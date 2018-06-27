from corna.constants import INTENSITY_COL, ISOTOPE_NA_MASS, KEY_NA
from corna.helpers import get_isotope_element, first_sub_second, parse_formula, chemformula_schema

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
    
    df['Label'] = [process_label(l) for l in df['Label']]
    
    s=df['Label'].apply(lambda x: x.split('_'))
    i=0
    for iso in iso_tracers:
        df[iso]= s.apply(lambda x: x[i+1])
        df[iso]= df[iso].astype(int)
        i=i+2

    del df['Label']
    return df

def process_corrected_df_for_metab(df, metabolite, isotracers, required_column):
    """
    This function returns formula dict for the formula of metabolite.
    formula- C5H4N3
    formula_dict- {'C':5,'H':4, 'N':3}
    """
    metabolite_df=df[df['Name']==metabolite]
    formula = metabolite_df.Formula.unique()
    formula_dict=parse_formula(formula[0])
    metabolite_df.set_index(isotracers, inplace=True) 
    required_df= metabolite_df.filter(required_column)
    return required_df, formula, formula_dict

def create_label_column(df, isotracers):
    """
    This function creates the label column back by joining information from multiple splited columns.
    Example:
            C13 N15 Label
            1   1   C13N15-Label-1-1 
    """
    df =df.rename_axis(None, axis=1).reset_index()
    isogroup=""
    df['Label']=""
    for i in isotracers:
        isogroup += i
        df['Label'] +='-'+ df[i].astype(str)
        del df[i]
    df['Label']= isogroup +'-label'+ df['Label']

    s=df['Label'].apply(lambda x: x.split('-'))

    return df

def add_info_to_final_df(info_df, metab, formula, iso_tracers):
    """
    Adds required columns back to the na corrected Dataframe 
    Required columns include : Label, Name, Formula
    """
    info_df=create_label_column(info_df, iso_tracers)
    info_df['Name']= metab
    info_df['Formula'] = formula
    return info_df

def save_original_df(df, iso_tracers):
    """
    This function modifies Label column from C13-Label-1 to C13N15-Label-1-0.
    And then returns back the original dataframe.
    """
    df= df.filter(['Name', 'Formula', 'Label', 'Sample', 'Intensity'])
    df['Original_label']=df['Label']
    df= convert_labels_to_std(df, iso_tracers)
    df= create_label_column(df, iso_tracers)
    df.drop('index', axis=1, inplace= True)
    return df
from model import Fragment
import helpers as hl
import numpy as np

def create_fragment_from_mass(name, formula, parent, isotope, isotope_mass, molecular_mass=None, mode=None):
    if molecular_mass != None:
        frag = Fragment(name, formula, parent, isotracer=isotope, isotope_mass=isotope_mass, molecular_mass=molecular_mass)
    elif mode != None:
        frag = Fragment(name, formula, parent, isotracer=isotope, isotope_mass=isotope_mass, mode=mode)
    else:
        frag = Fragment(name, formula, parent, isotracer=isotope, isotope_mass=isotope_mass)
    return {name:frag}

def create_fragment_from_number(name, formula, parent, label_dict):
    frag = Fragment(name, formula, parent, label_dict=label_dict)
    return {name:frag}

def create_combined_fragment(parent_fragment_dict, daughter_fragment_dict):
    parent_key, parent_fragment = parent_fragment_dict.ietms()[0]
    daughter_key, daughter_fragment = daughter_fragment_dict.items()[0]
    return {(parent_key, daughter_key): [parent_fragment, daughter_fragment]}

def validate_data(data):
    if not hl.check_if_all_elems_same_type(data.keys(), str):
        raise TypeError('Sample names should be of type string')
    if not hl.check_if_all_elems_same_type(data.values(), np.ndarray):
        raise TypeError('Intensities should be of type numpy arrays')

def add_data_fragment(fragment_dict, data):
    frag_key, frag = fragment_dict.items()[0]
    assert isinstance(data, dict)
    validate_data(data)
    return {frag_key : [frag, data]}

def add_data_isotopomers(frag_key, label_dict, intensity):
    label_key = hl.label_dict_to_key(label_dict)
    try:
        assert isinstance(intensity, np.ndarray)
    except AssertionError:
        raise AssertionError('intensity should be numpy array')
    label_intensity = {label_key: intensity}
    return {frag_key: label_intensity}

def parse_label_mass(label_mass):
    massdata = label_mass.split('_')
    try:
        isotracer = massdata[0]
        hl.get_isotope(isotracer)
        parent_mass = float(massdata[1])
        daughter_mass = float(massdata[2])
    except IndexError:
        raise IndexError('The key should have three components, isotope, parent mass '
                         'and daughter mass separated by _ in the same order')
    except KeyError:
        raise KeyError('First part of the key must be an isotope')
    except ValueError:
        raise ValueError('Masses should be convertible to floats')
    return {'tracer': isotracer, 'parent_mass': parent_mass, 'daughter_mass': daughter_mass}

def parse_label_number(label_number):
    return hl.create_dict_from_isotope_label_list(label_number.split('_'))

def parse_label_data_to_fragment(frag_data_dict, mass=False, number=False, mode=None):
    fragment_info, data = frag_data_dict.items()[0]
    label, sample_dict = data.items()[0]
    name = fragment_info[0]
    formula = fragment_info[1]
    parent = fragment_info[2]
    if mass == True:
            label_mass_dict = parse_label_mass(label)
            isotope = label_mass_dict['tracer']
            parent_mass = label_mass_dict['parent_mass']
            daughter_mass = label_mass_dict['daughter_mass']
            parent_frag = create_fragment_from_mass(name, formula, parent, isotope, parent_mass, mode=mode)
            daughter_frag = create_fragment_from_mass(name, formula, parent, isotope, daughter_mass, mode=mode)
            frag = create_combined_fragment(parent_frag, daughter_frag)
    elif number == True:
            label_number_dict = parse_label_number(label)
            frag = create_fragment_from_number(name, formula, parent, label_number_dict)
    else:
        raise TypeError('Labels in data should be deducable from mass or number')
    return add_data_fragment(frag, sample_dict)

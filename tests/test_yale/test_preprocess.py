import numpy
import pandas as pd

import corna.algorithms.background_correction as preproc
from corna.inputs.multiquant_parser import mq_df_to_fragmentdict
from corna.constants import ISOTOPE_NA_MASS

df=pd.DataFrame({'Sample':{0:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      1:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      2:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      3:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      4:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      5:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      6:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      7:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      8:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      9:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      10:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      11:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      },
                'Intensity':{0: 51666, 1:52361.856, 2:52538.436, 3:1874.364, 4:1483.272, 5:1874.364, 
                              6: 59689.272, 7:59950.872, 8:57204.072, 9:2746.8, 10:1525.128, 11:1438.8
                              },
                'Name':{0: 'dhap', 1:'dhap', 2:'dhap', 3:'dhap', 4:'dhap', 5:'dhap',
                        6: '2pg', 7:'2pg', 8:'2pg', 9:'2pg', 10:'2pg', 11:'2pg'
                        },
                'Component Name':{0: 'DHAP 169/97',1: 'DHAP 169/97', 2:'DHAP 169/97', 
                                  3: 'DHAP 170/97',4: 'DHAP 170/97',5: 'DHAP 170/97',
                                  6: '2PG 185/79',7: '2PG 185/79', 8: '2PG 185/79',
                                  9: '2PG 186/79', 10: '2PG 186/79', 11: '2PG 186/79'
                                  },
                'Formula':{0: 'H2O4P', 1: 'H2O4P', 2:'H2O4P', 3:'H2O4P',4:'H2O4P', 5:'H2O4P',
                           6: 'O3P',7: 'O3P',8: 'O3P',9: 'O3P',10: 'O3P',11: 'O3P'
                           },
                'Parent Formula':{0:'C3H6O6P', 1:'C3H6O6P',2: 'C3H6O6P',3: 'C3H6O6P',4:'C3H6O6P',5:'C3H6O6P',
                                  6:'C3H6O7P',7:'C3H6O7P',8:'C3H6O7P',9:'C3H6O7P',10:'C3H6O7P',11:'C3H6O7P'
                                  },
                'Label':{0:'C13_169.0_97.0',1:'C13_169.0_97.0',2:'C13_169.0_97.0',
                         3:'C13_170.0_97.0',4:'C13_170.0_97.0',5:'C13_170.0_97.0',
                         6:'C13_185.0_79.0',7:'C13_185.0_79.0',8:'C13_185.0_79.0',
                         9:'C13_186.0_79.0',10:'C13_186.0_79.0',11:'C13_186.0_79.0'
                         },
                'Background Sample':{0:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      1:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      2:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      3:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      4:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      5:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      6:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      7:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      8:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)',
                                      9:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                      10:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                      11:'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)'
                                      }                              

})

list_of_replicates = [numpy.array(['TA_SCS-ATP BCH_19May16_1June16.wiff (sample 1)',
                                           'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 17)',
                                           'TA_SCS-ATP BCH_19May16_1June16.wiff (sample 33)'
                                           ], dtype=object)]

background_df= pd.DataFrame({
    'Background Corrected':{0: 123.455, 1:-234.45}
})

def test_background_noise_label_daughter_unlabel():
    assert preproc.background_noise(31710, 0.011, 5, 1, 2, 0) == 1046.43

def test_background_noise_label_daughter_label():
     assert preproc.background_noise(31710, 0.011, 5, 1, 2, 1) == 697.62

def test_backround_subtraction():
    assert preproc.background_subtraction(27800, 2.29E+03) == 25510.0

def test_background_subtraction_negative():
     assert preproc.background_subtraction(278, 280) == 0

def test_background_correction():
     mq_fragdict = mq_df_to_fragmentdict(df)
     result = preproc.background_correction(mq_fragdict, list_of_replicates)
     output_list =[51666.0,52361.856 ,52538.436,
                           1720.4778000000001, 1329.3858, 1720.4778000000001,
                           59689.272, 59950.872,  57204.072,
                           1987.6527575999999, 765.9807575999996, 679.6527575999996
                           ]
    assert result['Background Corrected'].tolist() == output_list

def test_replace_negatives_background():
    result= preproc.replace_negatives_background(background_df)
    output_list= [123.455, 0]
    assert result['Background Corrected with zero'].tolist()== output_list

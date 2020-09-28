import pandas as pd

import corna.algorithms.nacorr_msms as algo
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
                'Intensity':{0: 51670, 1:52360, 2:52540, 3:1292.67, 4:901.67, 5:1292.67, 
                              6: 59690, 7:59950, 8:57200, 9:1969.77, 10:747.77, 11:661.77
                              },
                'Name':{0: 'DHAP 169/97', 1:'DHAP 169/97', 2:'DHAP 169/97', 3:'DHAP 169/97', 4:'DHAP 169/97', 5:'DHAP 169/97',
                        6: '2PG 185/79', 7:'2PG 185/79', 8:'2PG 185/79', 9:'2PG 185/79', 10:'2PG 185/79', 11:'2PG 185/79'
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
                         }                            

})

def test_na_correction_mimosa_without_background():
    test_df = algo.na_correction_mimosa(df, False)
   
    output_list= [53390.611000000004, -399.2437259999999, 61677.677, 
                25.82189399999993, 54103.588, -821.9009260000003,
                61946.33500000001, -1231.9645060000003, 54289.582, 
                -428.2147259999999, 59104.76000000001, -1228.2987060000003]
    
    assert test_df['NA Corrected'].tolist()== output_list
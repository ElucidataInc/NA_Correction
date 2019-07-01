'''
pyEQL Reaction Class

This file contains functions and methods for managing properties of 
individual chemical reactions in a solution.

Reaction objects are created on a per-Solution basis

:copyright: 2016 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

'''

## Dependencies
# import libraries for scientific functions
import math

# internal pyEQL imports
import pyEQL.activity_correction as ac
import pyEQL.water_properties as h2o
import pyEQL.solute as sol

# the pint unit registry
from pyEQL import unit

# import the parameters database
from pyEQL import paramsDB as db

# logging system
import logging
logger = logging.getLogger(__name__)

# add a filter to emit only unique log messages to the handler
import pyEQL.logging_system
unique = pyEQL.logging_system.Unique()
logger.addFilter(unique)

# add a handler for console output, since pyEQL is meant to be used interactively
ch = logging.StreamHandler()

# create formatter for the log
formatter = logging.Formatter('(%(name)s) - %(levelname)s - %(message)s')

# add formatter to the handler
ch.setFormatter(formatter)
logger.addHandler(ch)

class Reaction:
    '''
    Class representing a chemical reaction. Instances of this class 
    contain information about the products, reactants, equilibrium constant, etc.

    Parameters
    ----------
    reactants : dict
                Dictionary in which keys are strings representing the chemical formulas 
                of each reactant species and values are the corresponding stoichiometric coefficients.
    products : dict
                Dictionary in which keys are strings representing the chemical formulas 
                of each product species and values are the corresponding stoichiometric coefficients.
    log_K :         float
                Number representing the equilibrium constant of the reaction, dimensionless
    enthalpy:   Quantity, optional
                The enthalpy of the reaction, in units of energy/mol

    Returns
    -------
    Reaction
        A Reaction object.

    Examples
    --------

    See Also
    --------

    '''

    def __init__(self,reactants={},products={},log_K=0,**kwargs):
                            
        # set the products and reactants
        self.reactants=reactants
        self.products=products
        
        # set the equilibrium constant
        self.log_K = log_K
        
        # initialize optional keyword arguments
        self.delta_h = unit('0 J/mol')
        self.ref_temp = unit('25 degC')
        
        # if supplied, set the enthalpy
        if 'enthalpy' in kwargs:
            self.delta_h = unit(kwargs['enthalpy'])
        # if supplied, set the reference temperature for log K
        if 'ref_temp' in kwargs:
            self.ref_temp = unit(kwargs['ref_temp'])
        
    def Q(self,solution,activity_correction=True):
        '''
        Calculate the reaction quotient Q
        
        Parameters
        ----------
        activity_correction : bool, optional
                            If TRUE, the reaction quotient Q is calculated
                            using activities rather than concentraitons. If 
                            FALSE, mol/kg scale concentraitons are used.
                            Defaults to TRUE if omitted.
        '''
        numerator = denominator = 1
        
        #TODO what about pure liquids and solids?
        
        try:
            for species,coeff in self.products.items():
                if activity_correction is True:
                    numerator *= solution.get_activity(species) ** coeff
                else:
                    numerator *= solution.get_amount(species,'mol/kg').magnitude ** coeff
        except KeyError:
            # if one of the products isn't present, set it to zero
            numerator *= 0
        
        try:
            for species,coeff in self.reactants.items():
                if activity_correction is True:
                    denominator *= solution.get_activity(species) ** coeff
                else:
                    denominator *= solution.get_amount(species,'mol/kg').magnitude ** coeff
        except KeyError:
            # if one of the reactants isn't present, set it to zero
                denominator *= 0
        
        try:
            Q = numerator / denominator
        except ZeroDivisionError:
            # if the reactants aren't present, the reaction quotient is infinite
            import math
            Q = math.inf
            
        return Q
    
    def K(self,solution):
        '''
        Return the equilibrium constant for the reaction, adjusted for
        the solution temperature
        '''
        # if K was supplied as a parameter, use it
        K = 10 ** self.log_K
        
        # if K was not supplied, calculate a K from the Gibbs energy of formation
        # TODO
        
        # if the enthalpy is given and the solution temperature doesn't match the reference,
        # adjust the K for temperature
        if solution.get_temperature() != self.ref_temp:
            if self.delta_h.magnitude != 0:
                K_adjusted = self._adjust_temp_vanthoff(K,self.delta_h, solution.get_temperature(), self.ref_temp)
            else:
                logger.warning('Reaction enthalpy not provided, cannot adjust K for temperature.')
                K_adjusted = K
        else:
            K_adjusted = K
            
        return K_adjusted        
            
    def _error(self,solution):
        '''
        Helper method that calculates how far the reaction quotient Q is from K
        '''
        return abs(self.Q(solution)-self.K(solution.get_temperature()))
    
    def _adjust_temp_vanthoff(self,equilibrium_constant,enthalpy,temperature,reference_temperature):
        '''(float,float,number, optional number) -> float
        
        Helper method to adjust the reaction equilibrium constant from one temperature to another.
        
        Parameters
        ----------
        equilibrium_constant : float
                               The reaction equilibrium constant for the reaction
        enthalpy : Quantity
                   The enthalpy change (delta H) for the reaction in kJ/mol. Assumed
                   independent of temperature (see Notes).
        temperature : Quantity
                      the desired reaction temperature in degrees Celsius
        reference_temperature : Quantity
                          the temperature at which equilibrium_constant is valid.
       
        Returns
        -------
        float
            adjusted reaction equilibrium constant
        
        Notes
        -----
        This function implements the Van't Hoff equation to adjust measured 
        equilibrium constants to other temperatures. 
        
        .. math::
            ln(K2 / K1) = {\delta H \over R} ( {1 \over T_1} - {1 \over T_2} )
        
        This implementation assumes that the enthalpy is independent of temperature
        over the range of interest.[1]
        
        .. [1] Stumm, Werner and Morgan, James J. Aquatic Chemistry, 3rd ed, pp 53. 
            Wiley Interscience, 1996.
        
        Examples
        --------
        >>> adjust_temp_vanthoff(0.15,-197.6*unit('kJ/mol'),42*unit('degC'),25*unit('degC')) #doctest: +ELLIPSIS
        0.00203566...
        
        If the 'ref_temperature' parameter is omitted, a default of 25 C is used.
        
        >>> adjust_temp_vanthoff(0.15,-197.6*unit('kJ/mol'),42*unit('degC')) #doctest: +ELLIPSIS
        0.00203566...
        
        '''
        print(enthalpy)
        output = equilibrium_constant * math.exp( enthalpy / unit.R * ( 1 / reference_temperature.to('K') - 1 / temperature.to('K')))
        
        logger.info('Adjusted equilibrium constant K=%s from %s to %s degrees Celsius with Delta H = %s. Adjusted K = %s % equilibrium_constant,reference_temperature,temperature,enthalpy,output')
        
        logger.warning("Van't Hoff equation assumes enthalpy is independent of temperature over the range of interest")
        return output

    def __str__(self):
        #set output of the print() statement for the reaction
        str1 = ''
        counter=1
        length=len(self.reactants)
        for species,coeff in self.reactants.items(): 
            if coeff == 1:
                str1 += species
            else:
                str1 += str(coeff)+' '+species
            #only put pluses in between reactants, not after the last one
            if counter < length:
                str1 += ' + '
            counter += 1
            
        str1 += ' <=> '
        
        counter=1
        length=len(self.products)
        for species,coeff in self.products.items(): 
            if coeff == 1:
                str1 += species
            else:
                str1 += str(coeff)+' '+species
            # only put pluses in between reactants, not after the last one
            if counter < length:
                str1 += ' + '
            counter += 1

        return str1
    
    
    
    

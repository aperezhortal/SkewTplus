'''
New version of the thermodynamics.py module from the SkewT package:

https://pypi.python.org/pypi/SkewT

This module make use of the Cython functions in the _thermodynamics.pyx module
for Parcel analysis for computational efficiency efficiency 
and parallelization purposes
'''
# For python 3 portability
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from numpy import arctan, min, exp, log, where
import numpy
from numpy.ma.core import masked_array, masked_where, \
    masked_invalid, getmaskarray, getmask

from SkewTplus._thermodynamics import _parcelAnalysis4D, _parcelAnalysis3D, \
    _parcelAnalysis1D, _getLCL, _liftParcel, _moistAscent
from SkewTplus.errorHandling import fatalError


__author__ = "Andres Perez Hortal"
__copyright__ = "Copyright (c) 2017, Andres A. Perez Hortal"
__license__ = "BSD-3-Clause Licence, see LICENCE.txt for more details"
__email__ = "andresperezcba@gmail.com"
__credits__ = 'Based on the thermodynamics.py module from the SkewT package developed by Thomas Chubb'


#-----------------------------------------------------------------------
#   General Constants
#-----------------------------------------------------------------------
Rs_da = 287.05  # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v = 461.51  # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da = 1004.6  # Specific heat at constant pressure for dry air
Cv_da = 719.  # Specific heat at constant volume for dry air
Cp_v = 1870.  # Specific heat at constant pressure for water vapour
Cv_v = 1410.  # Specific heat at constant volume for water vapour
Cp_lw = 4218  # Specific heat at constant pressure for liquid water
Epsilon = 0.622  # Epsilon=Rs_da/Rs_v; The ratio of the gas constants
degCtoK = 273.15  # Temperature offset between K and C (deg C)
rho_w = 1000.  # Liquid Water density kg m^{-3}
grav = 9.80665  # Gravity, m s^{-2}
Lv = 2.5e6  # Latent Heat of vaporisation 
boltzmann = 5.67e-8  # Stefan-Boltzmann constant
mv = 18.0153e-3  # Mean molar mass of water vapor(kg/mol)
m_a = 28.9644e-3  # Mean molar mass of air(kg/mol)
Rstar_a = 8.31432  # Universal gas constant for air (N m /(mol K))
#-----------------------------------------------------------------------



def parcelAnalysis(pressure,
                   temperature,
                   dewPointTemperature,
                   hPa=False,
                   celsius=False,
                   maxPressureStep=2.,
                   fullFields=1,
                   method=0,
                   initialLevel=0,
                   useVirtual=1,
                   tolerance=0.1,
                   maxIterations=20,
                   depth=30000):
    """
    Function to perform a parcel analysis
    
    The shape of the input arrays (*pressure*, *temperature* and *dewPointTemperature*) 
    should be equal. 

    By default pressure and temperature are treated as Pascals and Kelvins degrees respectively.
    The working units can be selected by *hPa* and *celsius* keywords.
    The output values will automatically be converted to the selected units.
    
    The dimensions of the arrays should be one of the followings:
    
    * [time,height,latitude or y ,longitude or x] (3D + time)
    * [height,latitude or y ,longitude or x] (3D)
    * [height] 1D)
    
    See :py:func:`_parcelAnalysis1D` for a complete description of input parameters
    and returned values.  
    
    Parameters
    ----------
    
    hPa : bool, optional
        If True, pressure is treated as in hPa instead of Pa.
    
    celsius : bool, optional
        If True, temperature is considered as Celsius instead of Kelvin degrees.     
                          
    Other Parameters
    ----------------
    
    See :py:func:`_parcelAnalysis1D` for a complete description of input parameters.
    
    Returns
    -------
    
    See :py:func:`_parcelAnalysis1D` for a complete description of the returned values.
         
    """
    
    
    if ((pressure.shape != temperature.shape) or 
         (dewPointTemperature.shape != temperature.shape)):
        
        raise fatalError("Error computing Parcel Analysis",
                         "Input array shape mismatch",
                         "pressure.shape = %s" % str(pressure.shape),
                         "temperature.shape = %s" % str(temperature.shape),
                         "dewPointTemperature.shape = %s" % str(dewPointTemperature.shape),
                         )
    
    
    if isinstance(dewPointTemperature, masked_array):
        dewPointTemperature.data[dewPointTemperature.mask]=-9999.
        
    pressure = numpy.asarray(pressure, dtype=numpy.float32)
    temperature = numpy.asarray(temperature, dtype=numpy.float32)
    dewPointTemperature = numpy.asarray(dewPointTemperature, dtype=numpy.float32)
    
    if hPa:
        pressureConversionFactor = 100
        pressureUnits = 'hPa'
    else:
        pressureConversionFactor = 1 
        pressureUnits = 'Pa'
        
    if celsius:
        temperatureConversionOffset = degCtoK
        temperatureUnits = 'Celsius' 
    else:
        temperatureConversionOffset = 0
        temperatureUnits = 'Kelvin'
        
    args = (pressure*pressureConversionFactor, temperature+temperatureConversionOffset, dewPointTemperature+temperatureConversionOffset)
    kwargs = dict(maxPressureStep=maxPressureStep,
                  fullFields=fullFields,
                  method=method,
                  initialLevel=initialLevel,
                  tolerance=tolerance,
                  maxIterations=maxIterations,
                  useVirtual=useVirtual,
                  depth=depth)
                   
    if pressure.ndim == 4:
        myParcelAnalysis = _parcelAnalysis4D(*args, **kwargs)
    elif pressure.ndim == 3:
        myParcelAnalysis = _parcelAnalysis3D(*args, **kwargs)
    elif pressure.ndim == 1:
        myParcelAnalysis = _parcelAnalysis1D(*args, **kwargs)
    else:
        raise fatalError("Error computing Parcel Analysis",
                         "Number of array dimensions not supported"
                         "Number of array dimensions = %d" % pressure.ndim)
        
    
    
    # Convert the values to the selected units
    for key in myParcelAnalysis:
        
        if 'temperature' in key:
            myParcelAnalysis[key] -= temperatureConversionOffset
            
        if 'pressure' in key:
            myParcelAnalysis[key] /= pressureConversionFactor
            
        
    
    myParcelAnalysis['pressureUnits'] =   pressureUnits
    myParcelAnalysis['temperatureUnits'] =   temperatureUnits
      
    return myParcelAnalysis
        
    

def liftParcel(initialTemperature,
               pressureLevels,
               pressureAtLCL,
               initialLevel=0,
               hPa=False,
               celsius=False,
               maxPressureStep=1000):
    """
    Lift a parcel adiabatically.
    The parcel must correspond to a one of the pressure levels passed to the function. 
    
    By default pressure and temperature are treated as Pascals and Kelvins degrees respectively.
    The working units can be selected by *hPa* and *celsius* keywords.
    The temperature output values will automatically be converted to the selected units.
    
    .. _ndarray: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.htm
    
    
    Parameters
    ----------
    
    See :py:func:`SkewTplus._thermodynamics_liftParcel` for more  information on the parameters
    
    Returns
    -------
    
    parcelTemperature : ndarray_
        Parcel temperature at the pressure levels. All the levels below
        *initialLevel* are set to nan.
    
    """
    
    pressureLevels = numpy.asarray(pressureLevels, dtype=numpy.float32)
    pressureAtLCL = numpy.float32(pressureAtLCL)
    initialTemperature = numpy.float32(initialTemperature)
    
    if hPa:
        pressureConversionFactor = 100
    else:
        pressureConversionFactor = 1 
        
    
    if celsius:
        temperatureConversionOffset = degCtoK 
    else:
        temperatureConversionOffset = 0
        
    
    parcelTemperature = _liftParcel(initialTemperature + temperatureConversionOffset,
                                    pressureLevels * pressureConversionFactor,
                                    pressureAtLCL * pressureConversionFactor,
                                    initialLevel=initialLevel,
                                    maxPressureStep=1000)
    
    
    return parcelTemperature - temperatureConversionOffset     
    
    
    
def moistAscent(initialPressure,
                initialTemperature,
                finalPressure=-1,
                hPa=False,
                celsius=False,
                levels=500):

    """
    Lift a parcel moist adiabatically from initialPressure to finalPressure.
    By default pressure and temperature are treated as Pascals and Kelvins degrees respectively.
    The working units can be selected by *hPa* and *celsius* keywords.
    The output values will automatically be converted to the selected units.
    
    
    Parameters
    ----------
    
    initialPressure :  float
        Initial parcel pressure in Pascals by default
        
    initialTemperature :  float
        Initial temperature of the parcel in Kelvin degrees by default
        
    finalPressure : float32, optional
        Final parcel pressure in Pascals. Default 100hPa.
        
    hPa : bool, optional
        If True, pressure is treated as in hPa instead of Pa.
    
    celsius : bool, optional
        If True, temperature is considered as Celsius instead of Kelvin degrees.     
                          
    levels : int, optional
        Number of pressure levels between *initialPressure* and *finalPressure*.
        The pressure levels are spaced evenly on a log scale.
        
    Returns
    -------
    
    pressure : _ndarray, float32
        Pressure levels of the parcel
    
    temperature : _ndarray, float32
        Temperature of the parcel for the corresponding pressure levels
        
    """
    
    initialPressure = numpy.float32(initialPressure)
    initialTemperature = numpy.float32(initialTemperature)
    
    if hPa:
        pressureConversionFactor = 100
    else:
        pressureConversionFactor = 1 
        
    
    if finalPressure < 0:
        finalPressure = 100
        
    if celsius:
        temperatureConversionOffset = degCtoK 
    else:
        temperatureConversionOffset = 0
        
   
    pressure, temperature = _moistAscent(initialPressure=initialPressure * pressureConversionFactor,
                                         initialTemperature=initialTemperature + temperatureConversionOffset,
                                         finalPressure=finalPressure * pressureConversionFactor,
                                         levels=levels)
    
    return pressure / pressureConversionFactor, temperature - temperatureConversionOffset 
    
        
    
    
def getLCL(startPressure,
           startTemperature,
           startDewPointTemperature,
           hPa=False,
           celsius=False):
    """
    Get the Lifting Condensation Level of a parcel. 
    
    By default pressure and temperature are treated as Pascals and Kelvins degrees respectively.
    The working units can be selected by *hPa* and *celsius* keywords.
    The output values will automatically be converted to the selected units.  
    
    
    Parameters
    ----------
    
    startPressure :  float32
        Initial Pressure of the parcel in Pascals by default
    
        
    startTemperature :  float32
        Initial temperature of the parcel in Kelvin degrees by default
        
    
    startDewPointTemperature : float32
        Initial dew point temperature of the parcel in Kelvin degrees by default
        
    hPa : bool, optional
        If True, pressure is treated as in hPa instead of Pa.
    
    celsius : bool, optional
        If True, temperatures are considered as Celsius instead of Kelvin degrees.
        
        
    Returns
    -------
    
    pressure :  float32
        Pressure at the LCL 
    
    temperature : float32
        Temperature at the LCL         
        
    """
    
    
    
    startPressure = numpy.float32(startPressure)
    startTemperature = numpy.float32(startTemperature)
    
    if hPa:
        pressureConversionFactor = 100
    else:
        pressureConversionFactor = 1 
        
    
    if celsius:
        temperatureConversionOffset = degCtoK 
    else:
        temperatureConversionOffset = 0
        
        
    pressure, temperature = _getLCL(startPressure * pressureConversionFactor,
                                    startTemperature,
                                    startDewPointTemperature)
    
    return pressure / pressureConversionFactor, temperature - temperatureConversionOffset 

       
    
    
                     



def barometric_equation(presb_pa, tempb_k, deltah_m, Gamma=-0.0065):
    """The barometric equation models the change in pressure with 
    height in the atmosphere.

    INPUTS: 
    presb_k (pa):     The base pressure
    tempb_k (K):      The base temperature
    deltah_m (m):     The height differential between the base height and the 
                      desired height
    Gamma [=-0.0065]: The atmospheric lapse rate

    OUTPUTS
    pres (pa):        Pressure at the requested level

    REFERENCE:
    http://en.wikipedia.org/wiki/Barometric_formula
    """

    return presb_pa * (tempb_k / (tempb_k + Gamma * deltah_m)) ** (grav * m_a / (Rstar_a * Gamma))

def barometric_equation_inv(heightb_m, tempb_k, presb_pa, prest_pa, Gamma=-0.0065):
    """The barometric equation models the change in pressure with height in 
    the atmosphere. This function returns altitude given 
    initial pressure and base altitude, and pressure change.

    INPUTS: 
    heightb_m (m):
    presb_pa (pa):    The base pressure
    tempb_k (K)  :    The base temperature
    deltap_pa (m):    The pressure differential between the base height and the 
                      desired height

    Gamma [=-0.0065]: The atmospheric lapse rate

    OUTPUTS
    heightt_m

    REFERENCE:
    http://en.wikipedia.org/wiki/Barometric_formula
    """


    return heightb_m + tempb_k * ((presb_pa / prest_pa) ** (Rstar_a * Gamma / (grav * m_a)) - 1) / Gamma

def Theta(tempk, pres, pref=100000.):
    """Potential Temperature

    INPUTS: 
    tempk (K)
    pres (Pa)
    pref: Reference pressure (default 100000 Pa)

    OUTPUTS: Theta (K)

    Source: Wikipedia
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """

    try:
        minpres = min(pres)
    except TypeError:
        minpres = pres

    if minpres < 2000:
        print("WARNING: P<2000 Pa; did you input a value in hPa?")

    return tempk * (pref / pres) ** (Rs_da / Cp_da)

def TempK(theta, pres, pref=100000.):
    """Inverts Theta function."""

    try:
        minpres = min(pres)
    except TypeError:
        minpres = pres

    if minpres < 2000:
        print("WARNING: P<2000 Pa; did you input a value in hPa?")

    return theta * (pres / pref) ** (Rs_da / Cp_da)

def ThetaE(tempk, pres, e):
    """Calculate Equivalent Potential Temperature
        for lowest model level (or surface)

    INPUTS:
    tempk:      Temperature [K] 
    pres:       Pressure [Pa]
    e:          Water vapour partial pressure [Pa]

    OUTPUTS:
    theta_e:    equivalent potential temperature

    References:
    Eq. (9.40) from Holton (2004)
    Eq. (22) from Bolton (1980)
    Michael P. Byrne and Paul A. O'Gorman (2013), 'Land-Ocean Warming
    Contrast over a Wide Range of Climates: Convective Quasi-Equilibrium
    Theory and Idealized Simulations', J. Climate """

    # tempc
    tempc = tempk - degCtoK
    # Calculate theta
    theta = Theta(tempk, pres)


    # T_lcl formula needs RH
    es = VaporPressure(tempc)
    RH = 100.*e / es

    # theta_e needs q (water vapour mixing ratio)
    qv = vaporMixingRatio(e, pres)

    # Calculate the temp at the Lifting Condensation Level
    T_lcl = ((tempk - 55) * 2840 / (2840 - (log(RH / 100) * (tempk - 55)))) + 55

    # print "T_lcl :%.3f"%T_lcl

    #### DEBUG STUFF ####
    theta_l = tempk * (100000. / (pres - e)) ** (Rs_da / Cp_da) * (tempk / T_lcl) ** (0.28 * qv)
    # print "theta_L: %.3f"%theta_l

    # Calculate ThetaE
    theta_e = theta_l * exp((Lv * qv) / (Cp_da * T_lcl))

    return theta_e

def ThetaE_Bolton(tempk, pres, e, pref=100000.):
    """Theta_E following Bolton (1980)
    INPUTS:
    tempk:      Temperature [K] 
    pres:       Pressure [Pa]
    e:          Water vapour partial pressure [Pa]

    See http://en.wikipedia.org/wiki/Equivalent_potential_temperature
    """

    # Preliminary:
    T = tempk
    qv = vaporMixingRatio(e, pres)
    Td = DewPoint(e) + degCtoK
    kappa_d = Rs_da / Cp_da

    # Calculate TL (temp [K] at LCL):
    TL = 56 + ((Td - 56.) ** -1 + (log(T / Td) / 800.)) ** (-1)

    # print "TL: %.3f"%TL

    # Calculate Theta_L:
    thetaL = T * (pref / (pres - e)) ** kappa_d * (T / TL) ** (0.28 * qv)

    # print "theta_L: %.3f"%thetaL

    # put it all together to get ThetaE
    thetaE = thetaL * exp((3036. / TL - 0.78) * qv * (1 + 0.448 * qv))

    return thetaE

def ThetaV(tempk, pres, e):
    """Virtual Potential Temperature
    
    INPUTS
    tempk (K)
    pres (Pa)
    e: Water vapour pressure (Pa) (Optional)

    OUTPUTS
    theta_v    : Virtual potential temperature
    """ 

    mixr = vaporMixingRatio(e, pres)
    theta = Theta(tempk, pres)

    return theta * (1 + mixr / Epsilon) / (1 + mixr)

def GammaW(tempk, pres):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the environmental temperature and pressure.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Deg C/Pa)
    REFERENCE: 
    http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
    (Note that I multiply by 1/(grav*rho) to give MALR in deg/Pa)

    """

    tempc = tempk - degCtoK
    es = VaporPressure(tempc)
    ws = vaporMixingRatio(es, pres)

    # tempv=VirtualTempFromMixR(tempk,ws)
    tempv = VirtualTemp(tempk, pres, es)
    latent = Latentc(tempc)

    Rho = pres / (Rs_da * tempv)

    # This is the previous implementation:
    # A=1.0+latent*ws/(Rs_da*tempk)
    # B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    # Gamma=(A/B)/(Cp_da*Rho)

    # This is algebraically identical but a little clearer:
    A = -1.*(1.0 + latent * ws / (Rs_da * tempk))
    B = Rho * (Cp_da + Epsilon * latent * latent * ws / (Rs_da * tempk * tempk))
    Gamma = A / B

    return Gamma

def DensHumid(tempk, pres, e):
    """Density of moist air.
    This is a bit more explicit and less confusing than the method below.

    INPUTS:
    tempk: Temperature (K)
    pres: static pressure (Pa)
    mixr: mixing ratio (kg/kg)

    OUTPUTS: 
    rho_air (kg/m^3)

    SOURCE: http://en.wikipedia.org/wiki/Density_of_air
    """

    pres_da = pres - e
    rho_da = pres_da / (Rs_da * tempk)
    rho_wv = e / (Rs_v * tempk)

    return rho_da + rho_wv


def Density(tempk, pres, mixr):
    """Density of moist air

    INPUTS:
    tempk: Temperature (K)
    pres: static pressure (Pa)
    mixr: mixing ratio (kg/kg)

    OUTPUTS: 
    rho_air (kg/m^3)
    """
    
    virtualT = VirtualTempFromMixR(tempk, mixr)
    return pres / (Rs_da * virtualT)



def virtualTemp4(T, p):
    """Virtual Temperature, but using vapor mixing ratio instead the vapor pressure

    Parameters
    ---------    
    T: float
        Temperature in Celsius degrees
    
    Td: float
        Dew Point Temperature un Celsius Degrees
    
    p: float
        Atmospheric pressure (Pa)

    Returns
    -------
    Virtual temperature : float32
        In Celsius degrees
    """
    T = masked_invalid(T)
    _vaporMixingRatio = vaporMixingRatio(VaporPressure(T),p)
               
    return ((T+degCtoK) * (1. + (_vaporMixingRatio / Epsilon)) / (1. +_vaporMixingRatio)) - degCtoK

def virtualTemp3(T, Td, p):
    """Virtual Temperature, but using vapor mixing ratio instead the vapor pressure

    Parameters
    ---------    
    T: float
        Temperature in Celsius degrees
    
    Td: float
        Dew Point Temperature un Celsius Degrees
    
    p: float
        Atmospheric pressure (Pa)

    Returns
    -------
    Virtual temperature : float32
        In Celsius degrees
    """
    
    T = masked_invalid(T)
    Td = masked_where(numpy.array(Td)<-degCtoK,Td)
    Td = masked_invalid(Td)
    _vaporMixingRatio = vaporMixingRatio(VaporPressure(Td),p)
    
    T.filled(-99999)
    #print(_vaporMixingRatio)
    Tc = (T+degCtoK)
    factor = (1. + (_vaporMixingRatio.data / Epsilon)) / (1. +_vaporMixingRatio.data)
    
    _virtualTemp = (Tc*factor)  - degCtoK
    
    noTDMask = getmaskarray(_vaporMixingRatio)
    _virtualTemp.data[noTDMask]=T[noTDMask]
    _virtualTemp.mask = getmask(T)
    
    return _virtualTemp

def virtualTemp2(T, qv, p):
    """Virtual Temperature, but usig vapor mixing ratio instead the vapor pressure

    Parameters
    ---------    
    T: float
        Temperature in Kelvin degrees
    
    qv: float
        Vapour mixing ratio (Kg/Kg)
    
    p: float
        Atmospheric pressure (Pa)

    Returns
    -------
    Virtual temperature : float32
        In Kelvin degrees
    """

    return T * (1. + (qv / Epsilon)) / (1. + qv)  

def VirtualTemp(tempk, pres, e):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    e: vapour pressure (Pa)
    p: static pressure (Pa)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia)."""

    tempvk = tempk / (1 - (e / pres) * (1 - Epsilon))
    return tempvk
    

def VirtualTempFromMixR(tempk, mixr):
    """Virtual Temperature

    INPUTS:
    tempk: Temperature (K)
    mixr: Mixing Ratio (kg/kg)

    OUTPUTS:
    tempv: Virtual temperature (K)

    SOURCE: hmmmm (Wikipedia). This is an approximation
    based on a m
    """

    return tempk * (1.0 + 0.6 * mixr)

def Latentc(tempc):
    """Latent heat of condensation (vapourisation)

    INPUTS:
    tempc (C)

    OUTPUTS:
    L_w (J/kg)

    SOURCE:
    http://en.wikipedia.org/wiki/Latent_heat#Latent_heat_for_condensation_of_water
    """
   
    return 1000 * (2500.8 - 2.36 * tempc + 0.0016 * tempc ** 2 - 0.00006 * tempc ** 3)

def VaporPressure(tempc, phase="liquid"):
    """Water vapor pressure over liquid water or ice.

    INPUTS: 
    tempc: (C) OR dwpt (C), if SATURATION vapour pressure is desired.
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice

   
    RETURNS: e_sat  (Pa)
    
    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)
    
    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid = 6.112 * exp(17.67 * tempc / (tempc + 243.12)) * 100.
    over_ice = 6.112 * exp(22.46 * tempc / (tempc + 272.62)) * 100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase == "liquid":
        # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase == "ice":
        # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return where(tempc < 0, over_ice, over_liquid)
    else:
        raise NotImplementedError

def SatVap(dwpt, phase="liquid"):
    """This function is deprecated, return ouput from VaporPres"""

    print("WARNING: This function is deprecated, please use VaporPressure()" + 
            " instead, with dwpt as argument")
    return VaporPressure(dwpt, phase)



def vaporMixingRatio(e, p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure
          
    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon * e / (p - e)

def MixR2VaporPress(qv, p):
    """Return Vapor Pressure given Mixing Ratio and Pressure
    INPUTS
    qv (kg kg^-1) Water vapor mixing ratio`
    p (Pa) Ambient pressure
          
    RETURNS
    e (Pa) Water vapor pressure
    """

    return qv * p / (Epsilon + qv)


def DewPoint(e):
    """ Use Bolton's (1980, MWR, p1047) formulae to find tdew.
    INPUTS:
    e (Pa) Water Vapor Pressure
    OUTPUTS:
    Td (C) 
      """

    ln_ratio = log(e / 611.2)
    Td = ((17.67 - ln_ratio) * degCtoK + 243.5 * ln_ratio) / (17.67 - ln_ratio)
    return Td - degCtoK

def WetBulb(tempc, RH):
    """Stull (2011): Wet-Bulb Temperature from Relative Humidity and Air
    Temperature.
    INPUTS:
    tempc (C)
    RH (%)
    OUTPUTS:
    tempwb (C)
    """

    Tw = tempc * arctan(0.151977 * (RH + 8.313659) ** 0.5) + \
        arctan(tempc + RH) - arctan(RH - 1.676331) + \
        0.00391838 * RH ** 1.5 * arctan(0.023101 * RH) - \
        4.686035

    return Tw


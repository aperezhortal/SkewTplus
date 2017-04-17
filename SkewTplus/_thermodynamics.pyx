#!python
# cython: boundscheck=False
# cython: cdivision=True
# cython: embedsignature=True
'''
Module inspired in the thermodynamics.py from the SkewT package:

https://pypi.python.org/pypi/SkewT

The Parcel analysis was coded in Cython for computational efficiency efficiency
and parallelization purposes
'''

__author__ = "Andres Perez Hortal"
__copyright__ = "Copyright (c) 2017, Andres A. Perez Hortal"
__license__ = "BSD-3-Clause Licence, see LICENCE.txt for more details"
__email__ = "andresperezcba@gmail.com"


from cython.parallel import prange, parallel
from SkewTplus.errorHandling import fatalError

from cython.operator cimport dereference , preincrement
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf, fprintf, stderr

cimport cython

from libc.math cimport log, exp, ceil, log10

cimport numpy as numpy
import numpy as numpy

ctypedef numpy.float32_t float32
ctypedef numpy.float64_t float64


 
cdef float32 Epsilon = 0.622  # Epsilon=Rs_da/Rs_v; The ratio of the gas constants
cdef float32 degCtoK = 273.15  # Temperature offset between K and C (deg C)
cdef float32 T0 = 273.15  # Temperature offset between K and C (deg C)
cdef float32 Cp_da = 1004.6  # Specific heat at constant pressure for dry air
cdef float32 Rs_da = 287.05  # Specific gas const for dry air, J kg^{-1} K^{-1} 
cdef float32 Cp_ov_Rs_da = Cp_da / Rs_da 
cdef float32 Rs_ov_Cp_da = Rs_da / Cp_da
cdef float32 Rstar_a = 8.31432  # Universal gas constant for air (N m /(mol K))
cdef float32 g = 9.81
cdef float32 airMolarMass = 28.9644e-3  # Mean molar mass of air(kg/mol)


from common cimport float_abs, _linearInterpolation


cdef struct equilibriumPoint:
    float32 temperature  
    float32 pressure  
    float32 ABE
    int levelBelow
      
cdef struct temperatureAndPressure:
    float32 temperature  
    float32 pressure  
    int error    
    
cdef struct parcelAnalysisOutput:
    float32 pressureAtLCL
    float32 temperatureAtLCL
    float32 pressureAtLFC
    float32 temperatureAtLFC
    float32 pressureAtEL
    float32 temperatureAtEL
    float32 CIN 
    float32 CAPE
    int initialLevel
    int error 


# Python interface
def _moistLapseRateIntegrator(float32 initialTemperature,
                              float32 initialPressure,
                              float32 finalPressure,
                              float32 maxPressureStep=1000):
    """
    Integrate the moist lapse rate ODE , using a simple RK2 method, from an 
    initial pressure level to a final one. Only final step is returned.
    
    
    Parameters
    ----------
        
    initialTemperature :  float32
        Initial temperature in Kelvin degrees
        
    initialPressure :  float32
        Initial Pressure in Pa 
        
    finalPressure : float32
        Final Pressure in Pa.
        Temperature will be integrated up to this pressure level.
    
    maxPressureStep : float32, optional 
        Maximum pressure step allowed in the integration.
        If the difference between the initial and the final pressure levels
        exceeds this value, then the interval in divided in the necessary
        number of part in order to keep the pressure intervals below the threshold.    
        
    Returns
    -------
    temperature : float32
        Temperature at the final pressure level
        
    """
    return __moistLapseRateIntegrator(initialTemperature,
                                      initialPressure,
                                      finalPressure,
                                      maxPressureStep=maxPressureStep)
    


def _getLCL(float32 startPressure,
            float32 startTemperature,
            float32 startDewPointTemperature):
    """
    Get the Lifting Condensation Level (LCL) of a parcel
    
    Parameters
    ----------
    
    startPressure :  float32
        Initial Pressure of the parcel in Pascals
    
        
    startTemperature :  float32
        Initial temperature of the parcel in Kelvin degrees
        
    
    startDewPointTemperature : float32
        Initial dew point temperature of the parcel in Kelvin degrees
        
    Returns
    -------
    
    pressure :  float32
        Pressure at the LCL in Pa
    
    temperature : float32
        Temperature at the LCL in Kelvin degrees
    """
    
    cdef temperatureAndPressure myOutput  
    
    myOutput = _tempAndPressureAtLCL(startPressure,
                                     startTemperature,
                                     startDewPointTemperature)
    
    if myOutput.error > 0: 
        raise fatalError("Error computing the LCL")
     
    return myOutput.temperature , myOutput.pressure 

        
    
def _moistAscent(float32 initialPressure,
                 float32 initialTemperature,
                 float32 finalPressure=-1,
                 int levels=500):
    """
    Lift a parcel moist adiabatically from initialPressure to finalPressure.
    
    .. _ndarray: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.htm
    
    Parameters
    ----------
    
    initialPressure :  float32
        Initial parcel pressure in Pascals
        
    initialTemperature :  float32
        Initial temperature of the parcel in Kelvin degrees
        
    finalPressure : float32, optional
        Final parcel pressure in Pascals 
                          
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
    
    cdef float32 factor = 1.
    
    cdef numpy.ndarray[float32, ndim = 1] pressure
    cdef numpy.ndarray[float32, ndim = 1] temperature    
    
    temperature = numpy.zeros(levels, dtype=numpy.float32)
    
    
    if finalPressure < 0:
        finalPressure = 100 * 100  # 10 hPa
        
    pressure = numpy.logspace(log10(initialPressure), log10(finalPressure), levels, dtype=numpy.float32)
       
  
    temperature[0] = initialTemperature 
        
        
    cdef int i = 0
    
    for i in range(1, levels):
        
        temperature[i] = __moistLapseRateIntegrator(temperature[i - 1],
                                                    pressure[i - 1],
                                                    pressure[i])
        
    return pressure, temperature
    

def _liftParcel(float32 initialTemperature,
                numpy.ndarray[float32, ndim=1] pressure,
                float32 pressureAtLCL,
                int initialLevel=0,
                float32 maxPressureStep=1000):
    """
    Lift a parcel adiabatically.
    The parcel must correspond to a one of the pressure levels passed to the function. 
    
    
    .. _ndarray: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.htm
    
    
    Parameters
    ----------
        
    initialTemperature :  float32
        Initial temperature of the parcel in Kelvin degrees
        
    pressure :  ndarray_
        Pressure levels in Pascals where the parcel temperature is computed
    
    pressureAtLCL : float32
        Pressure at the LCL. 
        
                          
    initialLevel : int, optional
        Pressure level corresponding to the parcel.
        
    
    maxPressureStep : float32, optional
        Maximum pressure step (in Pascals) used in Moist Lapse rate integration.
        See :py:func:`_moistLapseRateIntegrator`.
     
        
    Returns
    -------
        
    parcelTemperature : ndarray_, float32
        Parcel temperature at the pressure levels. All the levels below
        *initialLevel* are set to nan.
    """   
    
    
    
    
    cdef int numberOfLevels = pressure.shape[0]
    cdef int level = 0
    cdef numpy.ndarray[float32, ndim = 1] parcelTemperature = numpy.zeros(numberOfLevels, dtype=numpy.float32)
    
    cdef float32 temperatureAtLCL = 0
    
    cdef float32  initialPressure = pressure[initialLevel]
    
    for level in range(initialLevel):
        parcelTemperature[level] = numpy.nan
    
    parcelTemperature[initialLevel] = initialTemperature
    
    for level in range(initialLevel, numberOfLevels - 1):

        if pressure[level + 1] >= pressureAtLCL:
            # Below LCL

            parcelTemperature[level + 1] = initialTemperature * ((pressure[level + 1] / initialPressure) ** Rs_ov_Cp_da)
            
        else:
            # Above LCL
            
            if pressure[level] > pressureAtLCL:

                temperatureAtLCL = initialTemperature * ((pressureAtLCL / initialPressure) ** Rs_ov_Cp_da)
                
                parcelTemperature[level + 1] = __moistLapseRateIntegrator(temperatureAtLCL,
                                                                       pressureAtLCL,
                                                                       pressure[level + 1],
                                                                       maxPressureStep=maxPressureStep)
                
            else:
                
                parcelTemperature[level + 1] = __moistLapseRateIntegrator(parcelTemperature[level],
                                                                       pressure[level],
                                                                       pressure[level + 1],
                                                                       maxPressureStep=maxPressureStep)
                
                
    return parcelTemperature     


        
def _parcelAnalysis4D(float32[:, :, :, :] pressure,
                      float32[:, :, :, :] temperature,
                      float32[:, :, :, :] dewPointTemperature,
                      float32 maxPressureStep=2.,
                      bint fullFields=True,
                      int method=0,
                      int initialLevel=0,
                      float32 tolerance=0.1,
                      int maxIterations=20,
                      float32 depth=30000):
    """
    Function to perform a parcel analysis over a 3D + time domain.
    
    The shape of the input arrays (*pressure*, *temperature* and *dewPointTemperature*) 
    should be equal and the dimensions must be
    ordered as [time,height,latitude or y ,longitude or x].
    
    Internally the arrays are treated as `Memory Views`_.
        
    .. _Memory Views: http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html
    
    See :py:func:`_parcelAnalysis1D` for a complete description of input parameters
    and returned values.  
         
    """ 
        
    cdef int i = 0, j = 0, k = 0, t = 0
    
    cdef int tSize = pressure.shape[0]  
    cdef int xSize = pressure.shape[3]  
    cdef int ySize = pressure.shape[2]
    cdef int zSize = pressure.shape[1]
           
      
    cdef float32[:, :, :] CIN
    cdef float32[:, :, :] CAPE

    CIN = numpy.zeros((tSize, ySize, xSize), dtype=numpy.float32)
    CAPE = numpy.zeros((tSize, ySize, xSize), dtype=numpy.float32)

            
    cdef float32[:, :, :] temperatureAtLCL, pressureAtLCL
    cdef float32[:, :, :] pressureAtLFC , temperatureAtLFC
    cdef float32[:, :, :] pressureAtEL, temperatureAtEL
    
        
    if fullFields:
        temperatureAtLCL = numpy.zeros_like(CIN)
        pressureAtLCL = numpy.zeros_like(CIN)
        
        temperatureAtLFC = numpy.zeros_like(CIN)
        pressureAtLFC = numpy.zeros_like(CIN)
        
        temperatureAtEL = numpy.zeros_like(CIN)
        pressureAtEL = numpy.zeros_like(CIN)
            
      
    cdef parcelAnalysisOutput myParcelAnalysis
    
    myParcelAnalysis.error = 0
         
    
    cdef int error = 0         
    
    for t in range(tSize):

        with nogil, parallel():
            
            for i in prange(xSize, schedule='dynamic'):
                
                for j in range(ySize):
                    
                    myParcelAnalysis = __parcelAnalysis1D(pressure[t, :, j, i],
                                                         temperature[t, :, j, i],
                                                         dewPointTemperature[t, :, j, i],
                                                         tolerance=tolerance,
                                                         maxIterations=maxIterations,
                                                         maxPressureStep=maxPressureStep,
                                                         method=method,
                                                         depth=depth,
                                                         initialLevel=initialLevel)
        
                    if myParcelAnalysis.error > 0:
                        with gil:
                            raise fatalError("Error during Parcel Analysis",
                                             "i=%d j%d" % (i, j))
                            
                             
                    CAPE[t, j, i] = myParcelAnalysis.CAPE
                    CIN[t, j, i] = myParcelAnalysis.CIN
                              
                    if fullFields:
                        temperatureAtLCL[t, j, i] = myParcelAnalysis.temperatureAtLCL
                        pressureAtLCL[t, j, i] = myParcelAnalysis.pressureAtLCL
                        
                        temperatureAtLFC[t, j, i] = myParcelAnalysis.temperatureAtLFC
                        pressureAtLFC[t, j, i] = myParcelAnalysis.pressureAtLFC
                        
                        temperatureAtEL[t, j, i] = myParcelAnalysis.temperatureAtEL
                        pressureAtEL[t, j, i] = myParcelAnalysis.pressureAtEL
    

            
    if fullFields:
        output = dict(temperatureAtLCL=numpy.asarray(temperatureAtLCL),
                      pressureAtLCL=numpy.asarray(pressureAtLCL),
                      temperatureAtLFC=numpy.asarray(temperatureAtLFC),
                      pressureAtLFC=numpy.asarray(pressureAtLFC),
                      temperatureAtEL=numpy.asarray(temperatureAtEL),
                      pressureAtEL=numpy.asarray(pressureAtEL),
                      initialLevel=initialLevel )
    else:
        output = dict()  
            
    output['CAPE'] = numpy.asarray(CAPE)
    output['CIN'] = numpy.asarray(CIN)  
                      
    return output    
        
   
def _parcelAnalysis3D(float32[:, :, :] pressure,
                      float32[:, :, :] temperature,
                      float32[:, :, :] dewPointTemperature,
                      float32 maxPressureStep=2.,
                      bint fullFields=True,
                      int method=False,
                      int initialLevel=0,
                      float32 tolerance=0.1,
                      int maxIterations=20,
                      float32 depth=30000):
    """
    Function to perform a parcel analysis over a 3D domain (at a single time)
    
    The shape of the input arrays (*pressure*, *temperature* and *dewPointTemperature*) 
    should be equal and the dimensions must be
    ordered as [height,latitude or y ,longitude or x].
    
    Internally the arrays are treated as `Memory Views`_.
        
    .. _Memory Views: http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html
    
    See :py:func:`_parcelAnalysis1D` for a complete description of input parameters
    and returned values.
    """
        
    cdef int i = 0, j = 0, k = 0
      
    cdef int xSize = pressure.shape[1]  
    cdef int ySize = pressure.shape[2]
    cdef int zSize = pressure.shape[0]
    
    cdef float32[:, :] CIN
    cdef float32[:, :] CAPE
            
    CIN = numpy.zeros((xSize, ySize), dtype=numpy.float32)
    CAPE = numpy.zeros((xSize, ySize), dtype=numpy.float32)
    
    cdef float32[:, :] temperatureAtLCL, pressureAtLCL
    cdef float32[:, :] pressureAtLFC , temperatureAtLFC
    cdef float32[:, :] pressureAtEL, temperatureAtEL
    
        
    if fullFields:
        temperatureAtLCL = numpy.zeros_like(CIN)
        pressureAtLCL = numpy.zeros_like(CIN)
        
        temperatureAtLFC = numpy.zeros_like(CIN)
        pressureAtLFC = numpy.zeros_like(CIN)
        
        temperatureAtEL = numpy.zeros_like(CIN)
        pressureAtEL = numpy.zeros_like(CIN)
            
      
    cdef parcelAnalysisOutput myParcelAnalysis
    
    myParcelAnalysis.error = 0
         
    
    cdef int error = 0         
    with nogil, parallel():
        for i in prange(xSize, schedule='dynamic'):
            
            if myParcelAnalysis.error > 0:                    
                break
            
            for j in range(ySize):
                
                myParcelAnalysis = __parcelAnalysis1D(pressure[:, i, j],
                                                     temperature[:, i, j],
                                                     dewPointTemperature[:, i, j],
                                                     maxPressureStep=maxPressureStep,
                                                     method=method,
                                                     initialLevel=initialLevel,
                                                     tolerance=tolerance,
                                                     maxIterations=maxIterations,
                                                     depth=depth)
                if myParcelAnalysis.error > 0:
                    with gil:
                        raise fatalError("Error during Parcel Analysis",
                                         "i=%d j%d" % (i, j))
                     
                         
                CAPE[i, j] = myParcelAnalysis.CAPE
                CIN[i, j] = myParcelAnalysis.CIN
                          
                if fullFields:
                    temperatureAtLCL[j, i] = myParcelAnalysis.temperatureAtLCL
                    pressureAtLCL[j, i] = myParcelAnalysis.pressureAtLCL
                    
                    temperatureAtLFC[j, i] = myParcelAnalysis.temperatureAtLFC
                    pressureAtLFC[j, i] = myParcelAnalysis.pressureAtLFC
                    
                    temperatureAtEL[j, i] = myParcelAnalysis.temperatureAtEL
                    pressureAtEL[j, i] = myParcelAnalysis.pressureAtEL
    

            
    if fullFields:
        output = dict(temperatureAtLCL=numpy.asarray(temperatureAtLCL),
                      pressureAtLCL=numpy.asarray(pressureAtLCL),
                      temperatureAtLFC=numpy.asarray(temperatureAtLFC),
                      pressureAtLFC=numpy.asarray(pressureAtLFC),
                      temperatureAtEL=numpy.asarray(temperatureAtEL),
                      pressureAtEL=numpy.asarray(pressureAtEL),
                      initialLevel=initialLevel 
                      )
    else:
        output = dict()  
            
    output['CAPE'] = CAPE
    output['CIN'] = CIN  
                      
    return output
      

def _parcelAnalysis1D(float32[:] pressure,
                      float32[:] temperature,
                      float32[:] dewPointTemperature,
                      float32 maxPressureStep=2.,
                      bint fullFields=True,
                      int method=0,
                      int initialLevel=0,
                      float32 tolerance=0.1,
                      int maxIterations=20,
                      float32 depth=30000):
    """
    Function to perform a parcel analysis for a vertical sounding (1D)
    
    The shape of the input arrays (*pressure*, *temperature* and *dewPointTemperature*) 
    should be equal and the dimension must correspond to height.
    
    Internally the arrays are treated as `Memory Views`_.
        
    .. _Memory Views: http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html
    
    The input arrays (*pressure*, *temperature* and *dewPointTemperature*) should be 
    one dimensional arrays.
        
    .. _Memory View: http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html
    
    .. _Memory Views: http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html
    
    
    Parameters
    ----------
        
    pressure :  `Memory View`_
        Pressure in Pascals
        
    temperature :  `Memory View`_
        Environmental temperature in Kelvin degrees
        
    dewPointTemperature :  `Memory View`_
        Environmental dew point temperature in Kelvin degrees
    

    maxPressureStep : float32, optional
        Maximum pressure step (in Pascals) used in Moist Lapse rate integration.
        See :py:func:`_moistLapseRateIntegrator`.
        

    method : int, optional
        Parcel analysis method used. Supported:
    
        * Most Unstable  : method=0
        * Single Parcel: method=1
        * Mixed Layer  : method=2 (Not supported yet)
             
         
    fullFields : bool, optional
        If False return only CAPE and CIN values. Otherwise temperature
        and pressure in the Lifting Condensation Level(LCL), the Level of
        Free Convection (LFC) and the Equilibrium Level (EL) is alse returned.
    
                          
    initialLevel : int, optional
        Initial level (index) used to compute the parcel analysis.
        Levels below this value are ignored.
        
        For the Single Parcel analysis, this level correspond to the parcel used.
        By default, the initial level is 0 (surface).
        
        For the Most Unstable method, this value is ignored.
        
    
    tolerance : float32 , optional
        Tolerance to be used in the cost function to find the LCL. 

    maxIterations : int, optional
        Maximum number of iterations to find the LCL using the bisection method.
        
        
    depth :  float32, optional
        Maximum depth in Pascals used to find the most unstable parcel.
         
        
        
    Returns
    -------
        
    parcelAnalysis : dict
        The parcel analysis return a dictionary with the results. Each element
        in the dictionary correspond to a float32 ndarray_ with the same shape
        as the input data arrays but with the vertical dimension removed.  
        
        The dictionary contains the following keys:
        
        * CAPE : Convective Available Potential Energy in Joules
         
        * CIN :  Convective Inhibition in Joules 
    
        if fullFields is True:
        
        * temperatureAtLCL :  Temperature at the LCL in Kelvins
        * pressureAtLCL : Pressure at the LCL in Pascals
        * temperatureAtLFC :  Temperature at the LFC in Kelvins
        * pressureAtLFC : Pressure at the LFC in Pascals
        * temperatureAtEL : Temperature at the EL in Kelvins
        * pressureAtEL : Pressure at the EL in Pascals

    """ 
    
    cdef parcelAnalysisOutput myParcelAnalysisOutput
    cdef int level = 0
    cdef float32 maxCAPE = 0 
    cdef float32 minpressure = pressure[initialLevel] - depth
    
    myParcelAnalysisOutput = __parcelAnalysis1D(pressure,
                                                temperature,
                                                dewPointTemperature,
                                                tolerance=tolerance,
                                                maxIterations=maxIterations,
                                                maxPressureStep=maxPressureStep,
                                                method=method,
                                                depth=depth,
                                                initialLevel=initialLevel)

    
    
        
        
    if fullFields:
        outputDict = dict(temperatureAtLCL=myParcelAnalysisOutput.temperatureAtLCL,
                          pressureAtLCL=myParcelAnalysisOutput.pressureAtLCL,
                          temperatureAtLFC=myParcelAnalysisOutput.temperatureAtLFC,
                          pressureAtLFC=myParcelAnalysisOutput.pressureAtLFC,
                          temperatureAtEL=myParcelAnalysisOutput.temperatureAtEL,
                          pressureAtEL=myParcelAnalysisOutput.pressureAtEL,
                          initialLevel=myParcelAnalysisOutput.initialLevel)
    else:
        outputDict = dict()  
        
    outputDict['CAPE'] = myParcelAnalysisOutput.CAPE
    outputDict['CIN'] = myParcelAnalysisOutput.CIN 
    
    return outputDict




cdef parcelAnalysisOutput __parcelAnalysis1D(float32[:] pressure,
                                             float32[:] temperature,
                                             float32[:] dewPointTemperature,
                                             float32 maxPressureStep=2.,
                                             int method=0,
                                             int initialLevel=0,
                                             float32 tolerance=0.1,
                                             int maxIterations=20,
                                             float32 depth=30000) nogil:                                            
    """
    Cython function used to perform a parcel analysis for a vertical sounding (1D).

    
    The input arrays (*pressure*, *temperature* and *dewPointTemperature*) should be 
    one dimensional arrays.
    
    This function do not use any type of Python objects to allow
    the release the GIL and called the function in for parallel processing 
    executions (http://cython.readthedocs.io/en/latest/src/userguide/external_C_code.html#releasing-the-gil).
    
    Internally the arrays are treated as Memory Views in order to allow the release of the GIL.
    See http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html for more information on memory views
    
    For a description of the Functions parameters see _parcelAnalysis1D function in this module.
     
    """  
    
    
     
    
    cdef parcelAnalysisOutput myParcelAnalysisOutput
    cdef parcelAnalysisOutput mostUnstableParcelAnalysis
    cdef int level = 0
    cdef float32 maxCAPE = 0 
    cdef float32 minpressure = 0
    cdef int MU_Level = 0 
    
    
    if method == 0:
        minpressure = pressure[0] - depth
        if minpressure < 0:
            minpressure = 0. 
        MU_Level = 0 
        myParcelAnalysisOutput = _singleParcelAnalysis1D(pressure,
                                                         temperature,
                                                         dewPointTemperature,
                                                         maxPressureStep=maxPressureStep,
                                                         tolerance=tolerance,
                                                         maxIterations=maxIterations,
                                                         initialLevel=0)
        
        # printf("Level:%d  CAPE:%.0f\n",level,myParcelAnalysisOutput.CAPE)
        level = 1   
        while pressure[level] > minpressure:
            
            mostUnstableParcelAnalysis = _singleParcelAnalysis1D(pressure,
                                                                 temperature,
                                                                 dewPointTemperature,
                                                                 maxPressureStep=maxPressureStep,
                                                                 tolerance=tolerance,
                                                                 maxIterations=maxIterations,
                                                                 initialLevel=level)
            
            # printf("Old CAPE=%.1f   CAPE=%.1f"%mostUnstableParcelAnalysis.CAPE , myParcelAnalysisOutput.CAPE)
            
            if  mostUnstableParcelAnalysis.CAPE > myParcelAnalysisOutput.CAPE:
                # printf("Level:%d  CAPE:%.0f\n",level,mostUnstableParcelAnalysis.CAPE)
                myParcelAnalysisOutput = mostUnstableParcelAnalysis
                # printf("Level:%d  CAPE:%.0f\n",level,myParcelAnalysisOutput.CAPE)
                MU_Level = level
                
            level += 1
            
        # printf("P(%d)=%.1f\n",MU_Level,pressure[MU_Level]/100)
    else:

        myParcelAnalysisOutput = _singleParcelAnalysis1D(pressure,
                                                         temperature,
                                                         dewPointTemperature,
                                                         tolerance=tolerance,
                                                         maxIterations=maxIterations,
                                                         maxPressureStep=maxPressureStep,
                                                         initialLevel=initialLevel)

    
    return  myParcelAnalysisOutput


                              

cdef parcelAnalysisOutput _singleParcelAnalysis1D(float32[:] pressure,
                                                  float32[:] temperature,
                                                  float32[:] dewPointTemperature,
                                                  float32 maxPressureStep=2.,
                                                  int maxIterations=20,
                                                  float32 tolerance=0.1,
                                                  int initialLevel=0) nogil:  
    """
    Cython function used to perform a single parcel analysis
    
    The input arrays (*pressure*, *temperature* and *dewPointTemperature*) should be 
    one dimensional arrays.
    
    For dewPointTemperature values, missing data is denoted by negative temperature values.
    
    This function do not use any type of Python objects to allow
    the release the GIL and called the function in for parallel processing 
    executions (http://cython.readthedocs.io/en/latest/src/userguide/external_C_code.html#releasing-the-gil).
    
    Internally the arrays are treated as Memory Views in order to allow the release of the GIL.
    See http://scipy-cookbook.readthedocs.io/items/ViewsVsCopies.html for more information on memory views
    
    For a description of the Functions parameters see _parcelAnalysis1D function in this module.
    
    IMPORTANT: The pressure levels are assumed decreasing in magnitude 
     
    """     
     
    cdef float32 temperatureAtLCL = 0 
    cdef float32 pressureAtLCL = 0 
    cdef parcelAnalysisOutput output  
    cdef float32 * parcelTemperature
    cdef float32 _auxPressure
    

    cdef int numberOfLevels = pressure.shape[0] 
    
    cdef int level = 0, next_level_LCL = numberOfLevels - 1  # next level after LCL
    cdef int lastLevel = numberOfLevels - 1  # Last level
    
    cdef float32 initialTemperature = temperature[initialLevel]
    cdef float32 initialPressure = pressure[initialLevel]
    cdef float32 initialDewPointTemperature = dewPointTemperature[initialLevel]
    
    if initialLevel >= numberOfLevels - 1:
        output.error = 2
        return output
        
    parcelTemperature = < float32 *> malloc(sizeof(float32) * numberOfLevels)
    for level in range(0, initialLevel):
        parcelTemperature[level] = 0 
    
    parcelTemperature[initialLevel] = initialTemperature
    
     
    #---------------------------------------------------------------------------
    # 1) First compute the parcel temperature up to the LCL ( dry ascent ) 
    #---------------------------------------------------------------------------
    
    # # Rough approximation
    cdef float32 J = 0
    cdef float32 initialMixingRatio = _vaporMixingRatio(_waterSatVaporPressureBolton(initialDewPointTemperature),
                                                        initialPressure)
    
    # TODO: check dew point temp is never higher than the env
    
    cdef bint useBisection = 0
    for level in range(initialLevel + 1, numberOfLevels):
        
        
        parcelTemperature[level] = initialTemperature * ((pressure[level] / initialPressure) ** Rs_ov_Cp_da)
        # printf("Level=%d   press=%.1f temp=%.1f \n",level,pressure[level],parcelTemperature[level] )        
        J = _T_LCL_RootFunction(parcelTemperature[level],
                                initialTemperature,
                                initialPressure,
                                initialMixingRatio
                                )
        
        # Lucky enough to find the root!
        if float_abs(J) < tolerance:
            # printf("Good , J=%.2f P(%d)=%.1f \n",J,level,pressure[level]/100)
            next_level_LCL = level + 1
            pressureAtLCL = pressure[level]
            temperatureAtLCL = parcelTemperature[level]
            useBisection = 0
            break
            
            
        if J > 0:  # Supersaturation
            next_level_LCL = level
            useBisection = 1
            break
        
        next_level_LCL = level + 1
        
        
    
    # # Use bisection method to improve the LCL value
    cdef int i = 0
    
    cdef float32 a 
    cdef float32 b 
    
   
    if useBisection:
        a = initialTemperature * ((pressure[next_level_LCL] / initialPressure) ** Rs_ov_Cp_da)
        J = _T_LCL_RootFunction(a,
                                initialTemperature,
                                initialPressure,
                                initialMixingRatio 
                                )
        
        b = parcelTemperature[next_level_LCL - 1]
        J = _T_LCL_RootFunction(b,
                                initialTemperature,
                                initialPressure,
                                initialMixingRatio 
                                )
    
    
        for i in range(maxIterations):
            
            temperatureAtLCL = (a + b) * 0.5
             
            J = _T_LCL_RootFunction(temperatureAtLCL,
                                    initialTemperature,
                                    initialPressure,
                                    initialMixingRatio 
                                    )
  
            if float_abs(J) < tolerance:
                break
            else:                
                if J < 0:
                    # Under-Saturated , T_LCL estimation -> b
                    b = temperatureAtLCL
                else:
                    a = temperatureAtLCL
    
    
    # printf("useBisection=%d  , J=%.2f   a=%.1f  b=%.1f\n",useBisection,J,a-273.15,b-273.15)                
    pressureAtLCL = initialPressure * ((temperatureAtLCL / initialTemperature) ** Cp_ov_Rs_da)     
            
    cdef float32 envTemperatureAtLCL        
    envTemperatureAtLCL = _linearInterpolation(pressureAtLCL,
                                               pressure[next_level_LCL - 1],
                                               pressure[next_level_LCL],
                                               temperature[next_level_LCL - 1],
                                               temperature[next_level_LCL])        
      
         
    #---------------------------------------------------------------------------                 
    # 2) Then compute moist ascent:
    # The moist ascent need to be computed by integrating a simple ODE
    #---------------------------------------------------------------------------
             
    
    if next_level_LCL < numberOfLevels:
        
        # Integrate from LCL to next level
        parcelTemperature[next_level_LCL] = __moistLapseRateIntegrator(temperatureAtLCL,
                                                                      pressureAtLCL,
                                                                      pressure[next_level_LCL])
        # Integrate for other levels
        if next_level_LCL + 1 < numberOfLevels:     
            for level in range(next_level_LCL + 1, numberOfLevels):
         
               
                parcelTemperature[level] = __moistLapseRateIntegrator(parcelTemperature[level - 1],
                                                            pressure[level - 1] ,
                                                            pressure[level])
                 
                # parcelTemperature[level] = lastTemperature
             
     
    #---------------------------------------------------------------------------
    # Now that we have the Temperature of the lifted parcel for each vertical level
    # 3) Compute Buoyancy (Using Virtual Temperature). B(P), not B(z)
    #---------------------------------------------------------------------------
     
    # accumulatedBuoyantEnergy calculated by a numerical integration
    # over the corresponding limits
    # Integrations performed using the trapezoidal rule 
    cdef float32 envVirtualTemp = 0 
    cdef float32 parcelVirtualTemp = 0
     
     
    cdef float32 envVaporMixingRatio = 0  
    cdef float32 parcelVaporMixingRatio = 0
     
    cdef float32 * buoyancy    
    buoyancy = < float32 *> malloc(sizeof(float32) * numberOfLevels)
    
    cdef float32 * accumulatedBuoyantEnergy    
    accumulatedBuoyantEnergy = < float32 *> malloc(sizeof(float32) * numberOfLevels)
    
    for level in range(0, initialLevel + 1):
        buoyancy[level] = 0 
        accumulatedBuoyantEnergy[level] = 0   
     
     
    # # First compute B between the surface and the LCL
    
    # # Vapor mixing ratio is constant 
    if next_level_LCL >= initialLevel + 1 :
         
        parcelVaporMixingRatio = _vaporMixingRatio(_waterSatVaporPressureBolton(dewPointTemperature[initialLevel]),
                                                   pressure[initialLevel])
        
        for level in range(initialLevel + 1, next_level_LCL):
            
            
            
            if dewPointTemperature[level] > 0: 
                envVirtualTemp = _virtualTemp(temperature[level],
                                              _waterSatVaporPressureBolton(dewPointTemperature[level]),
                                              pressure[level])
            else:
                envVirtualTemp = temperature[level]  # Missing data in humidity
                
            
            parcelVirtualTemp = _virtualTemp2(parcelTemperature[level], parcelVaporMixingRatio, pressure[level])
            
            # envVirtualTemp = temperature[level]
            # parcelVirtualTemp = parcelTemperature[level]            
             
            buoyancy[level] = Rs_da * (parcelVirtualTemp - envVirtualTemp) / pressure[level]
            
            accumulatedBuoyantEnergy[level] = (accumulatedBuoyantEnergy[level - 1] - 
                                              ( (pressure[level] - pressure[level - 1]) * 0.5 * 
                                                (buoyancy[level] + buoyancy[level - 1])) 
                                             ) 
             
     
    # # Then compute B above the LCL (saturation)    
    for level in range(next_level_LCL, numberOfLevels):
         
        if  dewPointTemperature[level] > 0:
            envVirtualTemp = _virtualTemp(temperature[level],
                                         _waterSatVaporPressureBolton(dewPointTemperature[level]) ,
                                         pressure[level])
        else:
            envVirtualTemp = temperature[level]
            
                                      
        parcelVirtualTemp = _virtualTemp(parcelTemperature[level],
                                         _waterSatVaporPressureBolton(parcelTemperature[level]),
                                         pressure[level])
        
        
         
        # envVirtualTemp = temperature[level]
        # parcelVirtualTemp = parcelTemperature[level]
        buoyancy[level] = Rs_da * (parcelVirtualTemp - envVirtualTemp) / pressure[level]
        
        accumulatedBuoyantEnergy[level] = (accumulatedBuoyantEnergy[level - 1] - 
                                           ((pressure[level] - pressure[level - 1]) * 0.5 * 
                                           (buoyancy[level] + buoyancy[level - 1])) 
                                           ) 
         
    
    #---------------------------------------------------------------------------   
    # 4) Find the Equilibrium levels where Buoyancy is zero and compute the
    # accumulated buoyant energy at each one of this levels
    #---------------------------------------------------------------------------
    

    # Find the LFC and the EL
    # EL: Equilibrium Level
    # LCL: Lifting condensation Level
    # LFC: Level of Free Convection
    # ABE Accumulated Buoyant Energy
    
    # The logic is as follows:
    # The EL correspond to the highest Equilibrium Point above the LCL with 
    # some positive buoyancy below it
    # Hence, we find the highest level with positive buoyancy above the LCL,
    # and then we find the EL between that level and the next
    # If the level is the Top, we set the Top as the ELs
    
    cdef int belowEL_level = -1
    cdef float32 pressureAtLFC = 0 , pressureAtEL = 0, accumulatedBuoyantEnergyAtLFC = 0
    cdef float32 temperatureAtLFC = 0 , temperatureAtEL = 0, accumulatedBuoyantEnergyAtEL = 0
    cdef float32  interpPressure = 0
    cdef vector[equilibriumPoint] _equilibriumPoints
    cdef bint EL_Found = 0  # True if an EL was found
    
    
    for level in range(next_level_LCL, numberOfLevels):
        if buoyancy[level] > 0:
            belowEL_level = level
    
    
    if belowEL_level >= 0:
        
        EL_Found = 1  
        
        if belowEL_level < lastLevel:        
            pressureAtEL = _linearInterpolation(0. ,
                                                buoyancy[belowEL_level], buoyancy[belowEL_level + 1],
                                                pressure[belowEL_level], pressure[belowEL_level + 1])
            
            temperatureAtEL = _linearInterpolation(pressureAtEL ,
                                                   pressure[belowEL_level], pressure[belowEL_level + 1],
                                                   parcelTemperature[belowEL_level], parcelTemperature[belowEL_level + 1])
            
            accumulatedBuoyantEnergyAtEL = _linearInterpolation(interpPressure ,
                                                              pressure[belowEL_level], pressure[belowEL_level + 1],
                                                              accumulatedBuoyantEnergy[belowEL_level], accumulatedBuoyantEnergy[belowEL_level + 1])
        
        else:
            pressureAtEL = pressure[lastLevel]
            temperatureAtEL = parcelTemperature[lastLevel]
            accumulatedBuoyantEnergyAtEL = accumulatedBuoyantEnergy[lastLevel]
            
                       
    
    
    # If no EL is found, pressure and temperature at EL and LFC are set to zero
    # (done already in variable initialization )
    
    cdef int minimumABE_Level = belowEL_level
    cdef float32 minimumABE = 1e10
    cdef equilibriumPoint myEquilibriumPoint
    cdef equilibriumPoint LFC_EquilibriumPoint
    
    cdef int point = 0 
    cdef vector[equilibriumPoint].iterator myIterator 
    
    if EL_Found:
        # If there is an EL, we define the LFC as the Equilibrium Point below the EL,
        # where accumulated buoyant energy is a minimum.
        
        # In an absolutely unstable state, the LFC could be below the LCL
        # This procedure differs from RIP CAPE one in:
        # the https://github.com/NCAR/wrf-python/blob/develop/fortran/rip_cape.f90 
        # because the LFC is allowed to be below the LCL
        
        # First we find the Equilibrium points
        for level in range(initialLevel + 1, belowEL_level):
         
            if float_abs(buoyancy[level]) < 0.00025 :
                myEquilibriumPoint.pressure = pressure[level]
                myEquilibriumPoint.temperature = parcelTemperature[level]
                myEquilibriumPoint.ABE = accumulatedBuoyantEnergy[level]
                
                _equilibriumPoints.push_back(myEquilibriumPoint)
                
            
            else:
                # If Buoyancy changes sign then there is an eq level between the levels
                if (buoyancy[level - 1] > 0) != (buoyancy[level] > 0):
                    
                    interpPressure = _linearInterpolation(0. ,
                                                          buoyancy[level - 1], buoyancy[level],
                                                          pressure[level - 1], pressure[level])
                    
                    myEquilibriumPoint.pressure = interpPressure
                
                    myEquilibriumPoint.temperature = _linearInterpolation(interpPressure ,
                                                                          pressure[level - 1], pressure[level],
                                                                          parcelTemperature[level - 1], parcelTemperature[level])
                    
                    myEquilibriumPoint.ABE = _linearInterpolation(interpPressure ,
                                                                  pressure[level - 1], pressure[level],
                                                                  accumulatedBuoyantEnergy[level - 1], accumulatedBuoyantEnergy[level])
            
            
                    _equilibriumPoints.push_back(myEquilibriumPoint)
                
        if _equilibriumPoints.size() > 0:
            # Finally we choose the Equilbrium point with the minimum ABE
            
            myIterator = _equilibriumPoints.begin()
            while myIterator != _equilibriumPoints.end():
                myEquilibriumPoint = dereference(myIterator)
                if myEquilibriumPoint.ABE < minimumABE:
                    LFC_EquilibriumPoint = myEquilibriumPoint
                    minimumABE = myEquilibriumPoint.ABE
            
                preincrement(myIterator)
            
            pressureAtLFC = LFC_EquilibriumPoint.pressure
            temperatureAtLFC = LFC_EquilibriumPoint.temperature
            accumulatedBuoyantEnergyAtLFC = LFC_EquilibriumPoint.ABE
        else:
            # If no EQ points below the EL are found but we have an EL, then ,
            # the parcel is absolutely unstable, and the LFC is equal to the
            # initial level.
            # There is no CIN present and Buoyancy is greater than zero in the
            # whole column (above the initial level )
            
            pressureAtLFC = pressure[initialLevel]
            temperatureAtLFC = temperature[initialLevel]
            accumulatedBuoyantEnergyAtLFC = 0                   
     
    #---------------------------------------------------------------------------
    # 5) Compute CAPE and CIN
    #---------------------------------------------------------------------------
              
              
    # # To compute CAPE and CIN perform a numerical integration 
    # # over the corresponding limits
    # # Integrations performed using the trapezoidal rule
    cdef float32 CAPE = 0.0 
    cdef float32 CIN = 0.0 
    
    if EL_Found:
        CAPE = accumulatedBuoyantEnergyAtEL - accumulatedBuoyantEnergyAtLFC
        CIN = -accumulatedBuoyantEnergyAtLFC
        
    
                     

    free(buoyancy) 
    free(parcelTemperature) 
    free(accumulatedBuoyantEnergy) 
    
    output.temperatureAtLCL = temperatureAtLCL
    output.pressureAtLCL = pressureAtLCL
    
    output.pressureAtLFC = pressureAtLFC
    output.temperatureAtLFC = temperatureAtLFC
    
    output.pressureAtEL = pressureAtEL
    output.temperatureAtEL = temperatureAtEL
    
    output.CAPE = CAPE
    output.CIN = CIN
    output.initialLevel = initialLevel
                    
    return output


     
      
cdef float32 __moistLapseRateIntegrator(float32 initialTemperature,
                                       float32 initialPressure,
                                       float32 finalPressure,
                                       float32 maxPressureStep=1000) nogil:
    """
    Integrate the moist lapse rate ODE , using a simple RK2 method, from an 
    initial pressure level to a final one. Only final step is returned.
    
    This function do not use any type of Python objects to allow
    the release the GIL and called the function in for parallel processing 
    executions (http://cython.readthedocs.io/en/latest/src/userguide/external_C_code.html#releasing-the-gil).
    
    See _moistLapseRateIntegrator function in this module for more details
    on the parameters and output.
    
        
    """
    
     
     
    cdef float32 dp, temperature , pressure , oldTemp , oldPres , k1
 
    cdef int numberOfInternalSteps, internalStep
 
     
    temperature = initialTemperature
    pressure = initialPressure
     
    dp = finalPressure - initialPressure
    numberOfInternalSteps = 1
     
     
    if float_abs(dp) > maxPressureStep:            
        numberOfInternalSteps = < int > ceil(float_abs(dp) / maxPressureStep)            
        dp = dp / (< float32 > numberOfInternalSteps) 
     
    for internalStep in range(numberOfInternalSteps):
         
        # Integrate from initialPressure to finalPressure
        # Runge - Kutta 2 
        k1 = dp * dTdP_Moist(temperature, pressure)
        temperature += dp * dTdP_Moist(temperature + k1 * 0.5, pressure + dp * 0.5)
        pressure += dp
                  
    return temperature 
 
 

     
  
cdef float32 dTdP_Moist(float32 temperature, float32 pressure) nogil:
    """
    Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the environmental temperature and pressure.
    
    The moist adiabatic lapse rate :math:`\Gamma = \frac{\partial T}{\partial P}` is computed as follows
     
    .. math::

       \Gamma = \frac{1+ \frac{L_v q_v}{R T }}{C_{pd} + \frac{L_v^2 q_v \epsilon}{R T^2}} \frac{R T_v}{P}  
    
    Where:

    * :math:`L_v` : latent heat of vaporization
    * :math:`C_{pd}` : specific heat at constant pressure of dry air
    * R : gas constant for dry air
    * T : temperature
    * :math:`T_v` : virtual temperature
    * :math:`\epsilon` : ratio of the gas constants for dry air and water vapor
    * :math:`q_v` : Vapor mixing ratio
        
    References:
        http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
          
 
    Parameters
    ---------  
    temperature : float32 
        Temperature in Kelvin degrees
    pressure : float 32 
        Pressure in Pascals (Pa)
     
 
    Returns
    -------
    
    Gamma: float32 
        The moist adiabatic lapse rate for the
        corresponding temperature and pressure (Deg K/Pa)
        
    """
     
    cdef float32 temperatureInC = temperature - degCtoK
     
    cdef float32 es = _waterSatVaporPressureBolton(temperature) 
    cdef float32 vaporMixingRatio = _vaporMixingRatio(es,
                                                      pressure)
     
    cdef float32 virtualTemp = temperature * (1. + (vaporMixingRatio / Epsilon)) / (1 + vaporMixingRatio)
     
     
    cdef float32 Lv = (2.501 - 0.00237 * temperatureInC) * 1e6  # Bolton, 1980 [J/kg]
     
     
    cdef float32 airDensity = pressure / (Rs_da * virtualTemp)
     
    cdef float32 A = 1.0 + Lv * vaporMixingRatio / (Rs_da * temperature) 
    cdef float32 B = Cp_da + ((Epsilon * Lv * Lv * vaporMixingRatio) / 
                                (Rs_da * temperature * temperature))
     
    cdef float32 Gamma = (A * Rs_da * virtualTemp) / (B * pressure) 
                            
 
    return Gamma


cdef inline float32 _T_LCL_RootFunction(float32 T ,
                                        float32 startTemperature,
                                        float32 startPressure,
                                        float32 startVaporMixingRatio) nogil:
    
    
    # Positive values means supersaturation  
    return (startVaporMixingRatio * startPressure * ((T / startTemperature) ** Cp_ov_Rs_da) / (startVaporMixingRatio + Epsilon) - 
            _waterSatVaporPressureBolton(T))
    



         
  
cdef temperatureAndPressure _tempAndPressureAtLCL(float32 startPressure,
                                                  float32 startTemperature,
                                                  float32 startDewPointTemperature) nogil:
     
    """
    Parameters
    ---------
     
    startPressure : float32
        Initial pressure in Pa
    startTemperature : float32
        Initial Temperature in Kelvin
         
    startDewPointTemperature : float32
        Initial Dew Point Temperature in Kelvin
    """
    
    cdef temperatureAndPressure myOutput  
    
    if float_abs(startTemperature - startDewPointTemperature) < 0.05:
        myOutput.temperature = startTemperature
        myOutput.pressure = startPressure        
        myOutput.error = 0 
        return myOutput
     
    cdef float32 startPressureInhPa = startPressure / 100
    cdef float32 startTemperatureInC = startTemperature - degCtoK
    cdef float32 startDewPointTemperatureInC = startDewPointTemperature - degCtoK
     
         
     
    cdef float32 startVaporMixingRatio = _vaporMixingRatio(_waterSatVaporPressureBolton(startDewPointTemperature),
                                                           startPressure)
     
    cdef float32 firstGuess = (startDewPointTemperature - 
                               (0.001296 * startDewPointTemperatureInC + 0.1963) * (startTemperature - startDewPointTemperature)
                               )
    
    # Bolton 1980 Approximation , use this value as an starting point to find the real LCL
    
    cdef float32 A = ((1. / (startDewPointTemperature - 56.)) + 
                       (log(startTemperature / startDewPointTemperature) / 800.))
                       
                        
    cdef float32 T_LCL = (1. / A) + 56.
    
    # Now use a bisection method to converge to the solution
    
      
     
    # cdef float32 T_LCL = newton(_T_LCL_RootFunction, firstGuess,
    #                            args=(startTemperature, startPressure, startVaporMixingRatio),
    #                            tol=1e-2)
     
    if T_LCL > startTemperature:
        printf("T_LCL=%f   , startTemperature = %f\n", T_LCL , startTemperature)
        myOutput.temperature = T_LCL
        myOutput.error = 1
        return myOutput
        
      
    cdef float32 pressure_LCL = startPressure * ((T_LCL / startTemperature) ** Cp_ov_Rs_da)
     
    myOutput.temperature = T_LCL
    myOutput.pressure = pressure_LCL    
    myOutput.error = 0     
    return myOutput
      
     
     
cdef inline float32 _waterSatVaporPressureBolton(float32 temperature) nogil:
    """
    Water vapor saturation pressure using Bolton 1980 formula.
    
    :math:`e_s = 6.112 exp^{\bigg( 17.67 T / (T + 243.5 )\bigg)}` 
    
    T : Temperature in Celsius 
    :math:`e_s` : Saturation vapor pressure in hPa
    
    Temperature is converted to Kelvins degrees internally
    
    Parameters
    ---------    
    temperature: float32 
        Temperature in Kelvin degrees 
        
    Returns
    -------
    
    saturationpressure : float32
        Water vapor saturation pressure in Pa
    
    """
    
    return 611.2 * exp(17.67 * (temperature - degCtoK) / (temperature + 243.5 - degCtoK))



cdef inline float32 _virtualTemp2(float32 T, float32 qv, float32 p) nogil:
    """Virtual Temperature, but usig vapor mixing ratio instead the vapor pressure

    Parameters
    ---------    
    T: float32 
        Temperature in Kelvin degrees
    
    qv: float32
        Vapour mixing ratio (Kg/Kg)
    
    p: float32
        Atmospheric pressure (Pa)

    Returns
    -------
    Virtual temperature : float32
        In Kelvin degrees
    """

    return T * (1. + (qv / Epsilon)) / (1. + qv)  




cdef inline float32 _virtualTemp(float32 T, float32 e, float32 p) nogil:
    """Virtual Temperature

    Parameters
    ---------    
    T: float32 
        Temperature in Kelvin degrees
    
    e: float32
        Vapour pressure (Pa)
    
    p: float32
        Atmospheric pressure (Pa)

    Returns
    -------
    Virtual temperature : float32
        In Kelvin degrees
    """

    return T / (1. - (e / p) * (1. - Epsilon))

    

cdef inline float32 _vaporMixingRatio(float32 e, float32 p) nogil:
    """
    Mixing ratio of water vapour. The input pressure must have the same units.
    
    Parameters
    ---------    
    e : float32 
        Vapor partial pressure , same units as total pressure
    
    p : total pressure 
        
    Returns
    -------
    
    vaporMixingRatio : float32
        Vapor Mixing ratio in kg/kg
        

    """

    return Epsilon * e / (p - e)



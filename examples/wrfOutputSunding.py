'''
Simple SkewT Plot of a Sounding from a WRF output file

For this example you need netCDF4

'''

# For python 3 portability
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from netCDF4 import Dataset
import numpy
from SkewTplus.sounding import sounding


#Load the WRF File
wrfOutputFile = Dataset("wrfOutputExample.nc")
theta = wrfOutputFile.variables["T"][:] + 300 # Potential temperature

# Pressure in hPa
pressure = (wrfOutputFile.variables['P'][:] + wrfOutputFile.variables['PB'][:]) 


qvapor = wrfOutputFile.variables['QVAPOR'][:]

qvapor = numpy.ma.masked_where(qvapor <0.00002, qvapor)

T0 = 273.15 
referencePressure = 100000.0  # [Pa]
epsilon = 0.622  # Rd / Rv

# Temperatures in Celsius
temperature = theta* numpy.power((pressure / referencePressure), 0.2854) - T0
vaporPressure = pressure * qvapor / (epsilon + qvapor)

dewPointTemperature = 243.5 / ((17.67 / numpy.log(vaporPressure * 0.01 / 6.112)) - 1.) #In celsius
dewPointTemperature = numpy.ma.masked_invalid(dewPointTemperature)


# Now we have the pressure, temperature and dew point temperature in the whole domain

# Select one vertical column , t =0 , x=30, y=30

inputData = dict(pressure=pressure[0,:,30,30]/100, 
                 temperature=temperature[0,:,30,30], 
                 dewPointTemperature=dewPointTemperature[0,:,30,30])

mySounding = sounding(inputData)
mySounding.quickPlot()
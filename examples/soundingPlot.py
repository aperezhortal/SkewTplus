'''
Simple SkewT Figure example from Sounding
'''

from __future__ import print_function, division

from SkewTplus.skewT import figure
from SkewTplus.sounding import sounding


#Load the sounding data
mySounding = sounding("./exampleSounding.txt")

# Create a Figure Manager 
mySkewT_Figure = figure()

# Add an Skew-T axes to the Figure
mySkewT_Axes = mySkewT_Figure.add_subplot(111, projection='skewx')

# Extract the data from the Sounding 
# The getCleanSounding method remove levels where invalid temperature or pressure
# values are present
pressure, temperature, dewPointTemperature = mySounding.getCleanSounding()

# Add a profile to the Skew-T diagram
mySkewT_Axes.addProfile(pressure,temperature, dewPointTemperature ,
                        hPa=True, celsius=True, method=0, diagnostics=True, useVirtual=1)
 
# Show the figure
mySkewT_Figure.show()


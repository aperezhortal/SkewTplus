'''
Simple SkewT Figure example from Sounding
'''

from __future__ import print_function, division


from SkewTplus.skewT import figure
from SkewTplus.sounding import sounding


#Load the sounding data
mySounding = sounding("./armSoundingExample.cdf")

# Create a Figure Manager 
mySkewT_Figure = figure()

# Add an Skew-T axes to the Figure
mySkewT_Axes = mySkewT_Figure.add_subplot(111, projection='skewx')

# Extract the data from the Sounding 
# The getCleanSounding method remove levels where invalid temperatura or pressure
# values are present
pressure, temperature, dewPointTemperature = mySounding.getCleanSounding()

# Add a profile to the Skew-T diagram
mySkewT_Axes.addProfile(pressure,temperature, dewPointTemperature ,
                        hPa=True, celsius=True, method=0, diagnostics=True)

mySkewT_Axes.set_title("El sondeo de la Paloma") 
# Show the figure
mySkewT_Figure.show()


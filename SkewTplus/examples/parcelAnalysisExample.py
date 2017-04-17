'''
Do a parcel analysis of the sounding and plot the parcel temperature
'''

from __future__ import print_function, division

from SkewTplus.skewT import figure
from SkewTplus.sounding import sounding
from SkewTplus.thermodynamics import parcelAnalysis, liftParcel


#Load the sounding data
mySounding = sounding("./exampleSounding.txt")

pressure, temperature, dewPointTemperature = mySounding.getCleanSounding()

# Perform a parcel analysis
# The full parcel analysis field is returned
# Most Unstable  parcel : method=0
# Start looking for the most unstable parcel from the first level (initialLevel=0)
# Use at maximum 5 iterations in the bisection method to find the LCL
# Since the sounding temperature and pressure are expressed in Celsius and hPa
# we set the corresponding keywords
myParcelAnalysis = parcelAnalysis(pressure,
                                  temperature,
                                  dewPointTemperature,
                                  hPa=True,
                                  celsius=True,
                                  fullFields=1,
                                  method=0,
                                  initialLevel=0,
                                  tolerance=0.1,
                                  maxIterations=20)

# Print the contents of the dictionary
for key,value in myParcelAnalysis.items():
    if isinstance(value, float) :
        print("%s = %.1f"%(key,value))
    else:
        print("%s = %s"%(key,str(value)))



#Plot the parcel trajectory in the SkewT diagram

# First we lift the parcel adiabatically
initialLevel = myParcelAnalysis['initialLevel']

parcelTemperature = liftParcel(temperature[initialLevel],
                               pressure,
                               myParcelAnalysis['pressureAtLCL'],
                               initialLevel=initialLevel,
                               hPa=True,
                               celsius=True)
                               
               
# Create a Figure Manager 
mySkewT_Figure = figure()

# Add an Skew-T axes to the Figure
mySkewT_Axes = mySkewT_Figure.add_subplot(111, projection='skewx')

# Plot the parcel temperature
mySkewT_Axes.plot(parcelTemperature, pressure, linewidth=3, color='r' )

# Add a marker for the LCL and the LFC
mySkewT_Axes.plot(myParcelAnalysis['temperatureAtLCL'], myParcelAnalysis['pressureAtLCL'], 
                  marker='o', color='b' , label='LCL')
mySkewT_Axes.plot(myParcelAnalysis['temperatureAtLFC'], myParcelAnalysis['pressureAtLFC'], 
                  marker='o', color='g' , label='LFC')

# Add a legend
mySkewT_Axes.legend(loc='center right')

mySkewT_Axes.set_title("Single Parcel Lifted adiabatically")

mySkewT_Figure.show()




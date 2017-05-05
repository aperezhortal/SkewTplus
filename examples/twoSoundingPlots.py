'''
Two soundings comparison, with the
parcel analysis, and plot them side to side
'''

from __future__ import print_function, division


from SkewTplus.skewT import figure
from SkewTplus.sounding import sounding



#Load the sounding data
mySounding1 = sounding("./bna_day1.txt")
mySounding2 = sounding("./bna_day2.txt")

# Create a Figure Manager with a suitable size for both plots
mySkewT_Figure = figure(figsize=(9,5))

# Now we want to create two axes side to side

# Add the first Skew-T axes to the Figure
mySkewT_Axes1 = mySkewT_Figure.add_subplot(121, projection='skewx',tmin=-40)


# Extract the data from the Sounding 
pressure, temperature, dewPointTemperature = mySounding1.getCleanSounding()


# Add a profile to the Skew-T diagram
mySkewT_Axes1.addProfile(pressure,temperature, dewPointTemperature ,
                        hPa=True, celsius=True, method=0, diagnostics=False)


mySkewT_Axes1.set_title("Day 1 Sounding")
 
# Add the second Skew-T axes to the Figure
mySkewT_Axes2 = mySkewT_Figure.add_subplot(122, projection='skewx',tmin=-40)

# Extract the data from the Sounding 
pressure, temperature, dewPointTemperature = mySounding2.getCleanSounding()

# Add a profile to the Skew-T diagram
mySkewT_Axes2.addProfile(pressure,temperature, dewPointTemperature ,
                        hPa=True, celsius=True, method=0, diagnostics=False)

mySkewT_Axes2.set_title("Day 2 Sounding") 

# Show the figure
mySkewT_Figure.show()


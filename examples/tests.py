'''
Tests
'''

# For python2.7 compatibility
from __future__ import print_function, division


from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy

from SkewTplus.skewT import figure
from SkewTplus.sounding import sounding
from SkewTplus.thermodynamics import parcelAnalysis, liftParcel
import matplotlib.pyplot as plt


# Read the examples soundings

# Load the sounding class from the package
#Load the sounding data

print("Testing TXT reader")
assert(sounding("./exampleSounding.txt")['pressure'][-1]==34.0)
assert(sounding("./94610.2010032200.txt")['pressure'][-1]==8.8) 
assert(sounding("./94866.2010030600.txt")['pressure'][-1]==37.6) 
assert(sounding("./94975.2013070900.txt")['pressure'][-1]==57.0)
assert(sounding("./94975.2013070900.txt")['pressure'][-1]==57.0)
assert(sounding("./bna_day1.txt")['pressure'][-1]==100.0)
assert(sounding("./bna_day2.txt")['pressure'][-1]==100.0)
assert(sounding("./sounding_high_tropo.txt")['pressure'][-1]==14.7)

# Show the SkewT diagram of the Sounding
mySounding = sounding("./exampleSounding.txt")
mySounding.quickPlot()


'''
Simple SkewT Figure example from Sounding
'''





#Load the sounding data
mySounding = sounding("./exampleSounding.txt")

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
 
# Show the figure
mySkewT_Figure.show()

'''
Plot two soundings in the same Skew-T diagram with out any parcel analysis
'''





#Load the sounding data
mySounding1 = sounding("./bna_day1.txt")
mySounding2 = sounding("./bna_day2.txt")

# Create a Figure Manager with a suitable size for both plots
mySkewT_Figure = figure(figsize=(5,6))

# Add the Skew-T axes to the Figure
mySkewT_Axes1 = mySkewT_Figure.add_subplot(111, projection='skewx',tmin=-40)


# Add one profile to the Skew-T diagram
# The line style is set to be a solid line and a label is added 
# to the plot. Since the label is not None, a legend will be added
# automatically to the plot
mySkewT_Axes1.addProfile(*mySounding1.getCleanSounding(),
                        hPa=True, celsius=True, parcel=False, 
                        label='Day 1', linestyle='-')


# Add a second profile to the Skew-T diagram
# The line style is set to be a dashed line 
# The location of the legend is specified to be 
# 'center right'
mySkewT_Axes1.addProfile(*mySounding2.getCleanSounding(),
                        hPa=True, celsius=True, parcel=False, 
                        label='Day 2', linestyle='--',loc='center right')

# Show the figure
mySkewT_Figure.show()

'''
Two soundings comparison, with the
parcel analysis, and plot them side to side
'''






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

'''
Do a parcel analysis of the sounding and plot the parcel temperature
'''




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



'''
Simple SkewT Figure example from Sounding


For this example you need netCDF4 and Basemap packages installed

To create the Data plots pyplot is used.
'''





#Load the WRF File
wrfOutputFile = Dataset("wrfOutputExample.nc")
theta = wrfOutputFile.variables["T"][:] + 300 # Potential temperature
pressure = wrfOutputFile.variables['P'][:] + wrfOutputFile.variables['PB'][:]


qvapor = wrfOutputFile.variables['QVAPOR'][:]

qvapor = numpy.ma.masked_where(qvapor <0.00002, qvapor)

T0 = 273.15 
referencePressure = 100000.0  # [Pa]
epsilon = 0.622  # Rd / Rv

temperature = theta* numpy.power((pressure / referencePressure), 0.2854) - T0

vaporPressure = pressure * qvapor / (epsilon + qvapor)

dewPointTemperature = 243.5 / ((17.67 / numpy.log(vaporPressure * 0.01 / 6.112)) - 1.) #In celsius
dewPointTemperature = numpy.ma.masked_invalid(dewPointTemperature)


# Now we have the pressure, temperature and dew point temperature in the whole domain
# Compute the parcel analysis for each vertical column and each time
#
# fullFields =0 , only return CAPE and CIN 
# Most Unstable  parcel : method=0
# Start looking for the most unstable parcel from the first level (initialLevel=0)
# Use at maximum 5 iterations in the bisection method to find the LCL
# Since the sounding temperature and pressure are expressed in Celsius and hPa
# we set the corresponding keywords
myParcelAnalysis = parcelAnalysis(pressure,
                                  temperature,
                                  dewPointTemperature,
                                  hPa=False,
                                  celsius=True,
                                  fullFields=0,
                                  method=0,
                                  initialLevel=0,
                                  tolerance=0.1,
                                  maxIterations=20)


## Create the Base Map for the CAPE color plot 

# Read the temperature and pressure fields
lon = wrfOutputFile.variables["XLONG"][0, :, :]
lat = wrfOutputFile.variables["XLAT"][0, :, :]


#---Read lat,lon for plotting
lon = wrfOutputFile.variables["XLONG"][0, :, :]
lat = wrfOutputFile.variables["XLAT"][0, :, :]


# Define and plot the meridians and parallels
min_lat = numpy.amin(lat)
max_lat = numpy.amax(lat)
min_lon = numpy.amin(lon)
max_lon = numpy.amax(lon)
    
# Create the basemap object
myBaseMap = Basemap(projection="merc",
                    llcrnrlat=min_lat,
                    urcrnrlat=max_lat,
                    llcrnrlon=min_lon,
                    urcrnrlon=max_lon,
                    resolution='h')
    
# Create the figure and add axes
myFigure = plt.figure(figsize=(8,8))
myAxes = myFigure.add_axes([0.1,0.1,0.8,0.8])

# Make only 5 parallels and meridians
parallel_spacing = (max_lat - min_lat) / 5.0
merid_spacing = (max_lon - min_lon) / 5.0
parallels = numpy.arange(min_lat, max_lat, parallel_spacing)
meridians = numpy.arange(min_lon, max_lon, merid_spacing)
    
myBaseMap.drawcoastlines(linewidth=1.5)
myBaseMap.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
myBaseMap.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

# Plot CAPE at time 0
CAPE = myParcelAnalysis['CAPE'][0,:]

myColorPlot = myBaseMap.pcolormesh(lon,lat, myParcelAnalysis['CAPE'][0,:],latlon=True, cmap='jet')
  
# Create the colorbar 
cb = myBaseMap.colorbar(myColorPlot,"bottom", size="5%", pad="5%")
cb.set_label("CAPE [J/kg]")
     
# Set the plot title
myAxes.set_title("CAPE")
           
plt.show()




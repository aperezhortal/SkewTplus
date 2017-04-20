.. _firstSteps:

Fist-Steps
==========


This chapter offers a quick overview of the main package capabilities:
To read a sounding from a txt file and create a quick SkewT diagram.


Sounding Data
-------------

The easiest way to get sounding data is to visit the University of 
Wyoming's website:

http://weather.uwyo.edu/upperair/sounding.html

To get some sounding data, follow the link and find the date, and location 
you are interested in, view the data as a text file and just save the file 
to your system. If you want to get loads of data please be considerate about 
the way you go about doing this! (Lots of wget requests makes the server 
administrators unhappy).

You can also pass your own data to SkewT by writing a text file in 
**identical** format to the University of Wyoming files, which are 
constant-width columns. Here's a sample of the first few lines of one of the 
bundled examples::

    94975 YMHB Hobart Airport Observations at 00Z 02 Jul 2013

    -----------------------------------------------------------------------------
       PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
        hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
    -----------------------------------------------------------------------------
     1004.0     27   12.0   10.2     89   7.84    330     14  284.8  306.7  286.2
     1000.0     56   12.4   10.3     87   7.92    325     16  285.6  307.8  286.9
      993.0    115   12.8    9.7     81   7.66    311     22  286.5  308.1  287.9


Alternatively you can create a dictionary with the column headers as keys 
and the data as 1D python arrays (preferably use ``ma.masked_array``). 
There's more about this under the "Running SkewT" section below.


Fist steps using SkewTplus
--------------------------

From now on, it's assumed that the package is installed and the current working
directory is the *examples* one, included in this package.

To read a sounding from a txt file and create a quick plot using the default
parameters we only have to do::

    from SkewTplus.sounding import Sounding
    
    #Load the sounding data
    mySounding = sounding("./exampleSounding.txt")
    mySounding.quickPlot()
        
The resulting plot will look like this:
        
.. image:: ../img/soundingQuickView.png

Now we can do the same thing, but with more control over the Figure:: 

    # Import the new figure class
    from SkewTplus.skewT import figure
    
    from SkewTplus.sounding import sounding
    
    
    #Load the sounding data
    mySounding = sounding("./exampleSounding.txt")
    
    # Create a Figure Manager 
    mySkewT_Figure = figure()
    
    # Add an Skew-T axes to the Figure
    mySkewT_Axes = mySkewT_Figure.add_subplot(111, projection='skewx')
    
    
    # Extract the data from the Sounding
    pressure = mySounding.soundingdata['pres']
    temperature =  mySounding.soundingdata['temp']
    dewPointTemperature = mySounding.soundingdata['dwpt']
    
    # Add a profile to the Skew-T diagram
    #  method=0 -> Most unstable parcel
    #  diagnostics -> add a text box in the Figure with the parcel analysis results
    mySkewT_Axes.addProfile(pressure,temperature, dewPointTemperature ,
                                hPa=True, celsius=True, method=0, diagnostics=True)
     
    # Show the figure
    mySkewT_Figure.show()

.. image:: ../img/soundingPlot.png

Lets now complicate the things a little bit and show one of the new capabilities
of the package. Let suppose that we want to compare two soundings, with the
parcel analysis, and plot them side to side::

    #Load the sounding data
    mySounding1 = sounding("./bna_day1.txt")
    mySounding2 = sounding("./bna_day2.txt")
    
    # Create a Figure Manager with a suitable size for both plots
    mySkewT_Figure = figure(figsize=(9,5))
    
    # Now we want to create two axes side to side
    
    # Add the first Skew-T axes to the Figure
    mySkewT_Axes1 = mySkewT_Figure.add_subplot(121, projection='skewx',tmin=-40)
    
    
    # Extract the data from the Sounding 
    pressure = mySounding1.soundingdata['pres']
    temperature =  mySounding1.soundingdata['temp']
    dewPointTemperature = mySounding1.soundingdata['dwpt']
    
    # Add a profile to the Skew-T diagram
    mySkewT_Axes1.addProfile(pressure,temperature, dewPointTemperature ,
                            hPa=True, celsius=True, method=0, diagnostics=False)
    
    
    mySkewT_Axes1.set_title("Day 1 Sounding")
     
    # Add the second Skew-T axes to the Figure
    mySkewT_Axes2 = mySkewT_Figure.add_subplot(122, projection='skewx',tmin=-40)
    
    # Extract the data from the Sounding 
    pressure = mySounding2.soundingdata['pres']
    temperature =  mySounding2.soundingdata['temp']
    dewPointTemperature = mySounding2.soundingdata['dwpt']
    
    # Add a profile to the Skew-T diagram
    mySkewT_Axes2.addProfile(pressure,temperature, dewPointTemperature ,
                            hPa=True, celsius=True, method=0, diagnostics=False)
    
    mySkewT_Axes2.set_title("Day 2 Sounding") 
    
    # Show the figure
    mySkewT_Figure.show()

.. image:: ../img/twoSoundingsPlots.png


The profile plotting capabilities are described in greater detail in the next chapter:
:ref:`profilePlotting`


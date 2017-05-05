'''
Created on Apr 22, 2017

@author: aperez
'''
from SkewTplus.sounding import sounding


mySounding = sounding()
mySounding.fetchFromWeb("20170410:00", "OAX")
mySounding.quickPlot()

mySounding = sounding("20170410:00",fileFormat='web', stationId= "OAX")
mySounding.quickPlot()


'''
Sounding class adapted from the SkewT python project Sounding,
developed by Thomas Chubb.

The UserDict base class was removed and changes for a dict one.

'''
from __future__ import division, print_function

from datetime import datetime

from numpy import ma
import numpy
from numpy.ma.core import mask_or, getmask, masked_invalid, masked_values

from SkewTplus.errorHandling import fatalError
from SkewTplus.skewT import figure


#TODO: Remove UserDict and used a dict class or change the class structure
class sounding(dict):
    # Copyright (c) 2013 Thomas Chubb 
    """Utilities to read, write and plot sounding data quickly and without fuss
    
    Parameters
    ----------
        
    filename:  str
        If creating a sounding from a file, the full file name. The 
        format of this file is quite pedantic and needs to conform 
        to the format given by the University of Wyoming soundings 
        (see weather.uwyo.edu/upperair/sounding.html) 
    data:  dict
        Soundings can be made from atmospheric data. This should be 
        in the form of a python dict with (at minimum) the following 
        fields (names case sensitive):

        temp: dry-bulb temperature (Deg C)
        dwpt: dew point temperature (Deg C)
        pres: pressure (hPa)
        sknt: wind speed (knots)
        wdir: wind direction (deg)

        The following fields are also used, but not required by the 
        plot_skewt routine:

        hght (m)
        relh (%)
        mixr (g/kg)
        thta (K)
        thte (K)
        thtv (K)
    """

    _AllowedKeys=['pres','hght','temp','dwpt','relh','mixr','drct','sknt','thta','thte','thtv']
    

    def __init__(self,filename=None,soundingdata=None):
        dict.__init__(self)

        self.soundingdata={}
        if soundingdata is None:
            self.uwyofile(filename)
        else:
            if not isinstance(soundingdata, dict):
                raise fatalError("Error initializing sounding instance",
                                 "soundingdata is not a dictionary",
                                 "type(soundingdata)=%s"%type(soundingdata))

            for key in soundingdata.keys():
                if key.lower() not in self._AllowedKeys:
                    self[key]=soundingdata.pop(key)
                else:
                    dd=soundingdata[key]
                    if hasattr(dd,'mask'):
                        ddm=dd
                    else:
                        ddm=ma.masked_invalid(dd)
                        ddm=ma.masked_values(ddm,-999)
                    ddm=ma.masked_array(ddm,mask=False).harden_mask()
                    self.soundingdata[key]=ddm
                if not self.has_key('StationNumber'): self['StationNumber']='(No Number)'
                if not self.has_key('SoundingDate'): self['SoundingDate']='(No Date)'

    
    
    def uwyofile(self,fname):
        """
        Reads the raw profile data from a University of Wyoming sounding file.

        This is the primary method of IO for SkewT. The University of 
        Wyoming maintains a nice database of global upper air data which is
        kept up-to-date. Given a filename, this method updates the sounding 
        data with the text data in the file.

        NOTES
        1. The input file has to conform *Exactly* to the University of 
           Wyoming file format. This is because I look for data fields at 
           specific places on each line.
        2. I ignore the diagnostics at the end of the file, because the idea 
           is to calculate these myself.
        3. When this no longer works I'll begin reading in a more array-esque 
           way.
        """
        #--------------------------------------------------------------------
        # This *should* be a convenient way to read a uwyo sounding
        #--------------------------------------------------------------------
        fid=open(fname)
        lines=fid.readlines()

        # New: handle whitespace at top of file if present
        while not lines[0].strip():
            lines.pop(0)


        nlines=len(lines)

        lhi=[1, 9,16,23,30,37,46,53,58,65,72]
        rhi=[7,14,21,28,35,42,49,56,63,70,77]

        # initialize output data structure
        output={}

        fields=lines[3].split()

        # Handle the file header
        # First line for WRF profiles differs from the UWYO soundings
        header=lines[0]
        if header[:5]=='00000':
            # WRF profile
            self['StationNumber']='00000'
            self['Longitude']=float(header.split()[5].strip(","))
            self['Latitude']=float(header.split()[6])
            self['SoundingDate']=header.split()[-1]
        else:
            self['StationNumber']=header[:5]
            dstr=(' ').join(header.split()[-4:])
            self['SoundingDate']=datetime.strptime(dstr,"%HZ %d %b %Y").strftime("%Y-%m-%d_%H:%M:%S") 

        # This is a data pre-initialisation step. I have used the number of lines minus the number
        # of lines of diagnostics.
        for ff in fields:
            # output[ff.lower()]=zeros((nlines-34))-999.
            output[ff.lower()]=[]

        lcounter=5
        for line,_ in zip(lines[6:],range(nlines)):
            lcounter+=1
            ### Version < 0.1.4
            # try: output[fields[0].lower()][idx]=float(line[lhi[0]:rhi[0]])
            # except ValueError: break

            ### New code. We test for pressure in the first column. 
            ### If there's no pressure, we get out!
            try: 
                output[fields[0].lower()].append(float(line[lhi[0]:rhi[0]]))
            except ValueError: 
                break
            
            for ii in range(1,len(rhi)):
                try: 
                    # Debug only:
                    # print fields[ii].lower(), float(line[lhi[ii]:rhi[ii]].strip())
                    ### Version < 0.1.4
                    # output[fields[ii].lower()][idx]=float(line[lhi[ii]:rhi[ii]].strip())
                     
                    ### New Code. Append to list instead of indexing 
                    ### pre-allocated data. Explicitly allocate -999
                    ### for invalid data (catch ValueError)
                    textdata=line[lhi[ii]:rhi[ii]].strip()
                    output[fields[ii].lower()].append(float(textdata))
                except ValueError: 
                    output[fields[ii].lower()].append(-9999.)

        for field in fields:
            ff=field.lower()
            # set mask for missing data
            
            dd=ma.masked_values(output[ff],-9999.)
            
            dd=ma.masked_array(dd,mask=False)
            
            dd.harden_mask()

            self.soundingdata[ff]=dd
            
            
        return None

    def getCleanSounding(self):
        """  
        Clean the sounding, when Temperature or pressure has a missing values,
        removing that pressure level """

        pressure = self.soundingdata['pres'].data
        temperature = self.soundingdata['temp'].data
        
        dewPointTemperature = self.soundingdata['dwpt']
        dewPointTemperature[dewPointTemperature.mask]=-9999
        dewPointTemperature = self.soundingdata['dwpt'].data
        
        mask = ~mask_or(getmask(self.soundingdata['pres']), getmask(self.soundingdata['temp']))
        
        return pressure[mask], temperature[mask], masked_values(dewPointTemperature[mask], -9999)
        
    
    def quickPlot(self, **kwargs):
        mySkewT_Figure = figure()

        # Add an Skew-T axes to the Figure
        mySkewT_Axes = mySkewT_Figure.add_subplot(111, projection='skewx',**kwargs)
        
        pressure, temperature, dewPointTemperature = self.getCleanSounding()
        
        # Add a profile to the Skew-T diagram
        mySkewT_Axes.addProfile(pressure, temperature, dewPointTemperature ,
                                hPa=True, celsius=True, method=0, diagnostics=True)
         
        # Show the figure
        mySkewT_Figure.show()
        
        



# For python 3 portability
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from builtins import (dict, object, range, str,
                      open, zip)
from datetime import datetime
import io
from os.path import isfile

from netCDF4 import Dataset
from numpy import ndarray
import numpy
from numpy.ma.core import mask_or, masked_values, \
    masked_invalid, getmaskarray, MaskedArray

from SkewTplus.errorHandling import fatalError
from SkewTplus.skewT import figure


# Long names from the TXT University of Wyoming soundings
_txtSoundingLongNames = {'pres' : 'Pressure',
                         'hght' : 'Height',
                         'alt' : 'Height',
                         'temp' : 'Temperature',
                         'dwpt' : 'Dew Point Temperature',
                         'relh' : 'Relative Humidity',
                         'mixr' : 'Vapor Mixing Ratio',
                         'drct' : 'Wind Direction',
                         'sknt' : 'Wind Speed',
                         'thta' : 'Potential Temperature',
                         'thte' : 'Equivalent Potential Temperature',
                         'thtv' : 'Virtual Potential Temperature'}


#: Dictionary with fields abbreviations ( abbreviation , full field name ) 
#: All the abbreviations are case insensitive.                   
abbreviations = {'pres' : 'pressure',
                 'hght' : 'height',
                 'alt' : 'height',
                 't' : 'temperature',
                 'temp' : 'temperature',
                 'dwpt' : 'dewPointTemperature',
                 'dp' : 'dewPointTemperature',
                 'relh' : 'relativeHumidity',
                 'rh' : 'relativeHumidity',
                 'mixr' : 'vaporMixingRatio',
                 'drct' : 'windDirection',
                 'sknt' : 'windSpeed',
                 'thta' : 'potentialTemperature',
                 'thte' : 'equivalentPotentialTemperature',
                 'thtv' : 'virtualPotentialTemperature',
                 'tdry' : 'temperature',
                 'wspd' : 'windSpeed',
                 'deg' : 'windDirection',
                 'v_wind' : 'vWind',
                 'u_wind' : 'uWind',
                 'asc' : 'ascentRate',
                 'lat' : 'latitude',
                 'lon': 'longitude',
                 'wstat' : 'windStatus' } 


    

class soundingArray(MaskedArray):
    """
    Class used to store the vertical sounding data arrays (MaskedArray_ with metadata)
    
    .. _MaskedArray: https://docs.scipy.org/doc/numpy/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray
    
    The class is a Masked array with extra attributes. These attributes are optional:
    
    For more information about masked arrays see:
    https://docs.scipy.org/doc/numpy/reference/maskedarray.html
    
    
    
    Attributes
    ----------
    
    longName : str
        Long name o the field value
    
    units : str    
        Units of the data values
        
    missingValue : float
        Missing data value. 
    """
         
    
    # Return a new instance of the class
    def __new__(cls, *args, **kwargs):
        # Input array is an already formed ndarray instance
                
                
        longName = kwargs.pop('longName', None)
        units = kwargs.pop('units', None)
        missingValue = kwargs.pop('missingValue', -9999.)
        
                
        obj = MaskedArray.__new__(cls, *args, **kwargs)
        
        # Add the new attributes to the created instance
        
        obj.longName = longName
        obj.units = units
        obj.missingValue = missingValue
        
        # Finally, we must return the newly created object:
        return obj
    





    
        

class sounding(object):
    """
    Class used to store Sounding data.
    
    It support quick initialization from txt files (WRF and University of Wyoming format),
    ARM sounding Netcdf files, or initialize from data in a dictionary.
    In addition the all the field can be loaded manually into the class by using the method
    :py:meth:`addField`.
    
    The data fields stored in the instance can be accessed like the elements of the dictionary
    using the [] (indexer) operators. For example::
    
        pressure = mySounding['pressure']
    
    In case that the field doesn't exist a :py:class:`~SkewTplus.errorHandling.fatalError`
    exception is raised.
    
    Several variable names abbreviations are supported. See :py:attr:`SkewTplus.sounding.abbreviations`.
    All the abbreviations are case insensitive.
    
    In addition, an iterator was implemented to access sequentially to all the fields
    and their corresponding values. For example to print all the fields stored and their corresponding 
    type is as simple as::
        mySounding = sounding('soundingExampleFile.txt')
        for fieldName, fieldValue in mySounding:
            print(fieldName , type(fieldValue))
    
    The list of the field names stored in the class can be retrieved with the 
    :py:meth:`fields` method.
    
    """
    
    missingValue = -9999.
    
    
    _soundingData = dict()
    
    def __init__(self, inputData=None, fileFormat=None):
        """
        Sounding Initialization
        
        It supports 3 types of initialization:
            
        * From a University of Wyoming txt file
        * Form an ARM sounding Netcdf file
        * From a dictionary with the field names, field values pairs
        
        All the field are added trough the :py:meth:`addField` method.
        The fields names support abbreviations when they are set or get.
        But internally the field names stored are not abbreviated.
        See :py:attr:`SkewT.sounding.abbreviations`.
        
        Parameters
        ----------
        
        inputData : str or dict
            If inputData is an string, it is considered a file path. 
            If the fileFormat keyword is not specified, the initialization will
            try to automatically detect the file format.
            If inputData is a dictionary, the sounding fields are taken from the
            dictionary (field names, field values pairs).
            If inputData is None, the create instance correspond to an empty sounding. 
            
        
        fileFormat : str
            File format. Supported values:
            * 'txt', "wrf", "uw or "wyoming" : University of Wyoming
            * "netcdf", "arm" : ARM netcdf soundings
            
            
        """
        
        self.addField('StationNumber', '(No Number)')
        self.addField('SoundingDate', '(No Date)')

        if inputData is not None:            
        
            if isinstance(inputData, str):
                # A file path was given
                
                
                if not isfile(inputData):
                    raise fatalError("The file does not exists: %s\n" % (inputData))
                
                if fileFormat is None:
                    # Try automatic detection of file format
                    
                    try:
                        self.readFromARMFile(inputData)
                    except OSError:
                        # If it is not a NETCDF , try TXT
                        self.readFromTxtFile(inputData)
                        
                else:   
                         
                    # If the argument is an string assume that is
                    # in the university of wyoming format 
                    if fileFormat.lower() in ('txt', "wrf", "uw", "wyoming") :
                        self.readFromTxtFile(inputData)
                    elif fileFormat.lower() in ('netcdf', 'arm'):
                        self.readFromARMFile(inputData)
                    else:
                        raise fatalError("Error loading Sounding",
                                         "Selected format not supported:%s" % fileFormat,
                                         "Accepted formats:%s" % str(('txt', "wrf", "uw", "wyoming", "netcdf")))
                        
            
            elif isinstance(inputData, dict):
                for fieldName, fieldValue in inputData.items():
                    
                    if isinstance(fieldValue, ndarray):
                        dataValues = masked_invalid(fieldValue)
                        dataValues = masked_values(fieldValue, self.missingValue)
                        dataValues.harden_mask()
                        self.addField(fieldName, fieldValue)
                    else:
                        self.addField(fieldName, fieldValue)
                                            
                
            else:
                raise fatalError("Input Data not supported",
                                 "type(inputData)=%s" % str(type(inputData)))
            
           
    def setMissingValue(self, missingValue):
        """ Set the missing data value"""
        self.missingValue = missingValue
    
    
        
    def _getFieldName(self, fieldName, exact=False):
        """
        Get the non-abbreviated file name. See 
        :py:attr:`abbreviations` for available abbreviations
        """
        
        if exact:
            fieldNameLowerCase = fieldName
        else:
            fieldNameLowerCase = fieldName.lower()
            
        if fieldNameLowerCase in abbreviations:
            fieldName = abbreviations[fieldNameLowerCase]
        
        return fieldName    
   
   
    
            
    def addField(self, fieldName, fieldValue, longName=None, units=None, **kwargs):
        """
        Add a field to the sounding
        
        .. _ndarray: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.htm
        
        Parameters
        ----------
        
        fieldName : str
            Name of the field to be added. It supports abbreviations 
            (see :py:attr:`abbreviations`)
        
        fieldValue : any
            Corresponding field value. If the field value is an ndarray_  then
            it will be converted to a :py:class:`~SkewTplus.sounding.soundingArray` with the corresponding 
            metadata.
            
        longName : str, optional
            Long name of the field
            
        units : str , optional
            units of the field
        
        missingValue : float, optional
            Value for missing data. Occurrences of this value will be masked.
            The default missingValue can be change with :py:meth:`setMissingValue` 
       
        """
        
        missingValue = kwargs.pop('missingValue', self.missingValue)
        
        fieldName = self._getFieldName(fieldName)   
        
        
        
        if isinstance(fieldValue, ndarray):
            
            _fieldValue = masked_invalid(fieldValue)
            _fieldValue = masked_values(_fieldValue , missingValue)
            _fieldValue.harden_mask()
            self._soundingData[fieldName] = soundingArray(_fieldValue,
                                                        longName=longName,
                                                        units=units,
                                                        missingValue=missingValue)
            
        else:
            self._soundingData[fieldName] = fieldValue
    
           
        
    def __getitem__(self, fieldName):
        """ Item getter """
        
        fieldName = self._getFieldName(fieldName)
        return self._soundingData[fieldName]
            
       
   
    def __iter__(self):
        """
        Return an iterator for the fields and corresponding values in the sounding
        """
        return iter(self._soundingData.items())
        
   
    
    def readFromTxtFile(self, filePath, headerLength=6):
        """
        Reads the raw profile data from a University of Wyoming sounding file.
        
        Notes
        -----
        1. Data can be downloaded from http://weather.uwyo.edu/upperair/sounding.html
        2. The input file has to conform *Exactly* to the University of 
           Wyoming file format. 
        3. The diagnostics at the end of the file are ignored.        
           
        """
        
        
        def isValidValue(valueString):
            # True if the string can be converted to float
            try: 
                _ = float(valueString)
                value = True
            except ValueError: 
                value = False
            return value
        
        
        with open(filePath, 'r') as fileId:
            
            lines = fileId.readlines()

            # New: handle whitespace at top of file if present
            offset=0
            while not lines[0].strip():
                lines.pop(0)
                offset +=1
            
            
            # Handle the file header
            # First line for WRF profiles differs from the UWYO soundings
            header = lines[0]
            if header[:5] == '00000':
                # WRF profile
                self.addField('StationNumber' , '00000')
                self.addField('Longitude' , float(header.split()[5].strip(",")))
                self.addField('Latitude' , float(header.split()[6]))
                self.addField('SoundingDate' , header.split()[-1])
            else:
                self.addField('StationNumber' , header[:5])
                dstr = (' ').join(header.split()[-4:])
                soundingDate = datetime.strptime(dstr, "%HZ %d %b %Y")
                self.addField('SoundingDate' , soundingDate.strftime("%Y-%m-%d_%H:%M:%S"))
            
            
            fieldsNames = lines[headerLength-3].split()
            filedsUnits = lines[headerLength-2].split()
            
            arePressureValues = [ isValidValue(_line[0:7]) for _line in lines]  
            
            fieldsData = numpy.genfromtxt(filePath, unpack=True,
                                          skip_header=headerLength+offset,
                                          max_rows=numpy.argwhere(arePressureValues).max()-headerLength+1,
                                          delimiter=7,
                                          usemask=True,
                                          missing_values = self.missingValue )
            

            for fieldName,fieldValues, units in zip(fieldsNames, fieldsData, filedsUnits):
                # Set mask for missing data
                                
                fieldValues = masked_invalid(fieldValues)                
                
                fieldValues.harden_mask()
                
                self.addField(fieldName , fieldValues,
                              units=units, missingValue=self.missingValue,
                              longName=_txtSoundingLongNames[fieldName.lower()])
            



    def readFromARMFile(self, filePath):
        """
        Reads the raw profile data from a ARM sounding file (Netcdf file).
        
        """
                
        armNetcdfFile = Dataset(filePath, mode="r")
        self.setMissingValue(-9999.)
        
        try:
            # The ARM field names must be in the _abbreviations attribute
            for field in ['pres', 'tdry', 'dp', 'wspd', 'deg', 'rh', 'u_wind',
                          'v_wind', 'wstat', 'asc', 'lat', 'lon', 'alt']:
                
                _variable = armNetcdfFile[field]
                if  hasattr(_variable, "missing_value"):
                    missingValue = _variable.missing_value
                else:
                    missingValue = -9999.
                
                self.addField(field , _variable[:],
                              units=_variable.units, missingValue=missingValue,
                              longName=_variable.long_name)
        except:
            raise
            
        finally:            
            armNetcdfFile.close()

        
    def iterFields(self):
        """ Return an iterator to all the fields stored"""
        return iter(self.__soundingData)
    
    
    def fields(self):
        """ Return a list with all the stored fields"""
        return self.__soundingData.keys()
    
        
    def getCleanSounding(self):
        """  
        Clean the sounding, when Temperature or pressure has a missing values,
        removing that pressure level """

        pressure = self['pressure']
        temperature = self['temperature']
        
        dewPointTemperature = self['dewPointTemperature']
        dewPointTemperature[getmaskarray(dewPointTemperature)] = self.missingValue
        
        mask = ~mask_or(getmaskarray(pressure), getmaskarray(temperature), shrink=False)

        return (pressure.data[mask],
                temperature.data[mask],
                masked_values(dewPointTemperature.data[mask], self.missingValue) 
                )
        
    
    
    
    def quickPlot(self, **kwargs):
        """ Do a quick plot of the sounding"""
        mySkewT_Figure = figure()

        # Add an Skew-T axes to the Figure
        mySkewT_Axes = mySkewT_Figure.add_subplot(111, projection='skewx', **kwargs)
        
        pressure, temperature, dewPointTemperature = self.getCleanSounding()
        
        # Add a profile to the Skew-T diagram
        mySkewT_Axes.addProfile(pressure, temperature, dewPointTemperature ,
                                hPa=True, celsius=True, method=0, diagnostics=True)
         
        # Show the figure
        mySkewT_Figure.show()
        
        


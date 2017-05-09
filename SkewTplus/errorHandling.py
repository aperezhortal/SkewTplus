"""
Error handling module
"""
# For python 3 portability
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

__author__ = "Andres Perez Hortal"
__copyright__ = "Copyright (c) 2017, Andres A. Perez Hortal, McGill University"
__license__ = "BSD-3-Clause License, see LICENCE.txt for more details"
__email__ = "andresperezcba@gmail.com"




class fileNotFoundException( Exception ):    
    """ Exception when a file is not found"""
    filePath = None
    message = None
    
    def __init__( self, filePath ):
        Exception.__init__( self, filePath ) 
        self.filePath = filePath
        self.message = 'Parameter file not found:  ' + self.filePath + '\n\n'

        
class fatalError( Exception ):   
    """ Fatal error exception """
       
    def __init__( self, mainErrorMessage ,*detailsMessages):
        """ Constructor """
        
        super(fatalError,self).__init__( self, mainErrorMessage ,*detailsMessages)
        
        self.mainErrorMessage = mainErrorMessage
        
        self.message = " " + mainErrorMessage + '\n\n'
        
        if len(detailsMessages) > 0:
            self.message += "Details:" + "\n"
        for message in detailsMessages:
            self.message += message + '\n'
        
        self.message +='\n'
            
    
    def __str__(self):
        return self.message    
    
    def __repr__(self):    
        return self.message    

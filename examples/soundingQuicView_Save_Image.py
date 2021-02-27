"""
Simple SkewT Figure example from Sounding
"""

# For python2.7 compatibility
from __future__ import print_function, division

# Load the sounding class from the package
from SkewTplus.sounding import sounding

# Load the sounding datas
mySounding = sounding("./exampleSounding.txt")

# If no X server is available to visualize the plot, save it to file
import matplotlib

matplotlib.use("Agg")

# Save the SkewT diagram of the Sounding
mySounding.quickPlot(output_file="quick_sounding.pdf")

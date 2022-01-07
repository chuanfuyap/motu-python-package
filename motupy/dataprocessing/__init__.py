greeting = "Hello World! This is a python package for processing and analysing microbiome data developed as part of the CURE project."

## OTUFile, the class object for reading in and processing OTU data
from .OTUdata import OTUdata
from .OTUnest import OTUnest

## GeneData, the class object for reading in and processing Gene data in microbiome currently only deal with humann2 data
from .GeneData import GeneData
from .GeneNest import GeneNest
greeting = "Hello world!"

## for testing of packagge installation
from .joke import joke

## OTUFile, the class object for reading in and processing OTU data
from .OTUdata import OTUdata
from .OTUnest import OTUnest

## describe
from .describe import describe

## methods to process reads dictionary, taken from OTUdata
from .mp_methods import *

from .kaiju_output import *
## OTUFile, the class object for reading in and processing OTU data
from motupy.dataprocessing.OTUdata import OTUdata
from motupy.dataprocessing.OTUnest import OTUnest

## GeneData, the class object for reading in and processing Gene data in microbiome currently only deal with humann2 data
from motupy.dataprocessing.GeneData import GeneData
from motupy.dataprocessing.GeneNest import GeneNest

__all__ = ["OTUdata", "OTUnest", "GeneData","GeneNest"]
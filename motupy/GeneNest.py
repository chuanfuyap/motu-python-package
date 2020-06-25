"""
    GeneNest object is for storing and processing multiple sample of gene data.
    The data will be stored as a dictionary {key: sample_ID, value: GeneData object}
"""
import pandas as pd 
from .GeneData import GeneData
import os 

class GeneNest:
    def __init__(self, verbose = 0):
        """
        Creates GeneNest object for manipulation/transformation

        Parameters
        ------------
        N/A

        Returns
        ------------
        N/A

        """ 
        self.verbose = verbose
        self.nest = {}

    def build_from_folder(self, input_folder, input_type, extension=None):
        """
        Creates OTUnest object for manipulation/transformation

        Parameters
        ------------
        input_folder: str,
            folder location of the input files (taxa profiling output file)

        input_type: str,
            what type of taxa profiling tool was used to generate the file, "E.g. Kaiju"

        Returns
        ------------
        N/A

        """ 
        self.nest = {}

        files = [f for f in os.listdir(input_folder) if extension in f]
        files_dir = [input_folder+'/'+f for f in files]
        for f, n in zip(files_dir, files):
            if self.verbose>=1 : print(n)
            tmpGeneData = GeneData(f, input_type, self.verbose, extension)
            
            self.nest[tmpGeneData.file_id] = tmpGeneData.genefile

        return self.to_dataframe()

    def to_dataframe(self):
        """
        Turns OTUnest nest object from dictionary to a pandas dataframe.

        Parameters
        ------------
        N/A 

        Returns
        ------------
        pandas dataframe object

        """ 
        return pd.DataFrame(self.nest).transpose()

"""
    OTUnest object is for storing and processing multiple sample of OTU data.
    The data will be stored as a dictionary {key: sample_ID, value: OTUData object}
"""
import pandas as pd 
from OTUdata import OTUdata
import os 
class OTUnest:
    def __init__(self, verbose = 0):
        """
        Creates OTUnest object for manipulation/transformation

        Parameters
        ------------
        N/A

        Returns
        ------------
        N/A

        """ 
        self.verbose = verbose
        self.basic_ranks = ['superkingdom',
                        'phylum',
                        'class',
                        'order',
                        'family',
                        'genus',
                        'species']
        self.ranks = {'superkingdom': set(),
                        'phylum': set(),
                        'class': set(),
                        'order': set(),
                        'family': set(),
                        'genus': set(),
                        'species': set()}
        self.superkingdom = {'virus' : set(),
                            'bacteria' : set(),
                            'eukaryote' : set(),
                            'archaea' : set()}

        
        self.nest = {}

    def build_from_folder(self, input_folder, input_type, extension=None, artifact_threshold=0, make_clade_relative=True, cumulate=False):
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
            tmpOTUdata = OTUdata(f, input_type, artifact_threshold, self.verbose, extension)
            
            if make_clade_relative:
                tmpOTUdata.turn_reads_to_clade_relative_abundance()
            elif cumulate:
                tmpOTUdata.taxa_cumulation()
                
            for rank in self.ranks.keys():
                self.ranks[rank].update(tmpOTUdata.ranks[rank])
            for sk in self.superkingdom.keys():
                self.superkingdom[sk].update(tmpOTUdata.superkingdom[sk])

            self.nest[tmpOTUdata.file_id] = tmpOTUdata.otufile

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

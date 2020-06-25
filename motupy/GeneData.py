"""
    GeneData object is for storing and processing single sample of gene data, be it metabolic potential (DNA) or metabolic activity (RNA)
    The data will be stored as a dictionary {key: Gene ID, value: read counts/relative abundance}
"""
from .humann2_output import *
from .describe import describe

class GeneData:
    def __init__(self, file_loc, input_type, verbose=0, extension=None):
        """
        Creates OTUdata object for manipulation/transformation

        Parameters
        ------------
        file_loc: str,
            location/file name of the input file (taxa profiling output file)

        input_type: str,
            what type of taxa profiling tool was used to generate the file, "E.g. Humann2"

        Returns
        ------------
        N/A

        """ 

        self.file_id = ''
        self.verbose = verbose

        if extension is None:
            self.file_id = self.process_file_name(file_loc)
        else:
            self.file_id = self.process_file_name_with_known_extension(file_loc, extension)

        if input_type.lower()  not in ['humann2']:
            print("[ERROR] only 'humann2_mergedgenefamilies input type supported")
        if input_type.lower() == 'humann2':
            self.genefile = read_h2(file_loc)


    def process_file_name(self, file_loc):
        """
        Parameters
        ------------
        file_loc : str,
            the name/location of the file being processed

        Returns
        ------------
        id_name: str,
            the id_name to be associated with GeneData object
        """
        slash_split = file_loc.split('/')
        slash_split = slash_split[-1]
        dot_split = slash_split.split('.')
        
        id_name = dot_split[0]
        
        return id_name

    def process_file_name_with_known_extension(self, file_loc, ext):
        """
        Parameters
        ------------
        file_loc : str,
            the name/location of the file being processed
        
        ext : str,
            the known/desired extension used to cut off the name for the id generation

        Returns
        ------------
        id_name: str,
            the id_name to be associated with GeneData object
        """
        slash_split = file_loc.split('/')
        slash_split = slash_split[-1]
        id_name = slash_split.replace(ext, '')
        
        return id_name
"""
    This module contains methods for Exploratory Data Analysis (EDA). 
    
    As this is built on 'classic' python methods for EDA, it accepts DataFrame object where row=sample, column=features (in our case OTU).
"""

import numpy as np
import pandas as pd

from ete3 import NCBITaxa; ncbi = NCBITaxa()

from scipy.spatial import distance
from skbio.stats import composition, ordination
from skbio import DistanceMatrix

class EDA():
    def describe_df(self, dataframe):
        """
        TO DO
        Parameters
        ------------
        dataframe : pandas DataFrame, 
            the desired taxanomy_ranks made up of str

        ranks_dict : dict,
            dictionary containing:
                key = taxanomic ranks
                value = taxa id belonging to the given rank

        Returns
        ------------
        N/A
        """
        ## this should print out summary statistics of the dataframe, 
        ## sample total, OTU total,read average for each taxonomic rank

        df = dataframe.copy()
        ranks = self.update_ranks(df)
        
        print('Sample Count\t\t|%d'%df.shape[0])
        print('Total OTU Count\t\t|%d'%df.shape[1])
        
        stats_df = pd.DataFrame()

        for key in ranks.keys():
            stats = df[ranks[key]].sum(0).describe()
            stats_df[key.capitalize()] = stats

        return stats_df


    def filter_otu_coverage(self, dataframe, otu_coverage=0.0025):
        """
        Selects OTU that are found in given fraction of samples. 
        
        Parameters
        ------------
        dataframe : pandas DataFrame,
                        rows = sample,
                        columns = OTU

        otu_coverage : float,
                        the fraction of samples covered by the the OTU

        Returns
        ------------
        dataframe : pandas DataFrame that has been filtered
        """
        df = dataframe.copy()
        df.replace(0, np.nan, inplace=True)
        df.dropna(axis=1, thresh=(otu_coverage*df.shape[0]), inplace =True)
        df.fillna(0, inplace=True)

        return df

    def clean_artifact(self, dataframe, artifact_threshold=5):
        """
        Replaces count values less than artifact_threshold values with zeroes.

        Idea is that values that found to be 1 could be an artifact from kaiju. 
        
        Parameters
        ------------
        dataframe : pandas DataFrame,
                        rows = sample,
                        columns = OTU

        artifact_threshold : int,
                        the threshold cut off for potential artifacts.

        Returns
        ------------
        dataframe : pandas DataFrame that has been updated 
        """
        pass

    def update_ranks(self, dataframe):
        """
        This method collects all the taxa id in the dataframe's columns and group them into their respective taxanomy rank 
        and store them in a dictionary (ranks). 

        Parameters
        ------------
        dataframe: pandas dataframe object
            where 
                rows = samples
                columns = OTU
                and each datapoint is the either read count or relative abundance. 

        Returns
        ------------
        ranks: dict object
            where {key: taxanomy rank, value: a set containing all the taxanomy id beloninging to the rank}
        """
        df = dataframe.copy()
        ranks = {'superkingdom': set(),
                            'phylum': set(),
                            'class': set(),
                            'order': set(),
                            'family': set(),
                            'genus': set(),
                            'species': set()}
        for taxid in df.columns:
            taxid = int(taxid)
            try:
                taxid_rank = ncbi.get_rank([taxid])[taxid] 
            except:
                print(taxid,'has no rank')
            if taxid_rank in ranks:
                ranks[taxid_rank].add(str(taxid))
        return ranks

    def get_taxid_sets(self, dataframe):
        """
        This method collects all the taxa id in the dataframe's columns and group them into their respective taxanomy rank 
        and store them in a dictionary separate set objects. 

        Parameters
        ------------
        dataframe: pandas dataframe object
            where 
                rows = samples
                columns = OTU
                and each datapoint is the either read count or relative abundance. 

        Returns
        ------------
        species, genus, family, order, class, phylum, superkingdom: set objects 
            they all contain tax id belonging to each taxonomic rank.
        """
        df = dataframe.copy()
        ranks = {'superkingdom': set(),
                            'phylum': set(),
                            'class': set(),
                            'order': set(),
                            'family': set(),
                            'genus': set(),
                            'species': set()}
        for taxid in df.columns:
            taxid = int(taxid)
            try:
                taxid_rank = ncbi.get_rank([taxid])[taxid] 
            except:
                print(taxid,'has no rank')
            if taxid_rank in ranks:
                ranks[taxid_rank].add(str(taxid))

        return ranks['species'], ranks['genus'], ranks['family'], ranks['order'], ranks['class'], ranks['phylum'], ranks['superkingdom']
        
    def aitchison_distance_matrix(self, df):
        """
        This function takes the a count dataframe data and does the following:
            i) impute the zeroes with composition.multiplicative_replacement.
            ii) centred log-ratio transformation on data
            iii) calculate the aitchison distance (which is the euclidean distance on clr transformed values)
        
        Parameters
        ------------
        dataframe: pandas dataframe object
            where 
                rows = samples
                columns = OTU
                and each datapoint is in read count 

        Returns
        ------------
        dataframe: pandas dataframe object
            this is a distance matrix stored as in the pandas dataframe,
            where rows=columns, the samples, with a diagonal of zeroes.
        """
        copied = df.copy()

        copied.fillna(0, inplace=True)
        X_imputed = composition.multiplicative_replacement(copied)
        X_clr = composition.clr(X_imputed)
        
        ##note the dataframe is transposed here, so the columns are now samples, 
        # making the follow functions iterate over the column for subject comparisons
        copied = pd.DataFrame(X_clr, columns=copied.columns, index=copied.index).T

        col_names = copied.columns
        
        scored = {}
        
        ## fill up dictionary with empty dicts for each column names
        for col_1 in col_names:
            scored[col_1] = {}
        
        for col_1 in col_names:
            for col_2 in col_names:
                distance_score = distance.euclidean(copied[col_1], copied[col_2])
                scored = self.add_to_dict(scored, col_1, col_2,  distance_score)
                scored = self.add_to_dict(scored, col_2, col_1,  distance_score)
                
        return pd.DataFrame(scored)

    def add_to_dict(self, the_dict, column_name, index_name, value):
        """
        Called during the generation of distance matrix to store the pairwise sample distance scores to prevent repeats from being calculated and stored. 
        
        Parameters
        ------------
        the_dict: dict object
        column_name: str
        index_name: str
        value: float

        Returns
        ------------
        the_dict: dict object
            contains all the calculated value so far. 

        """
        if index_name not in the_dict[column_name]:
            the_dict[column_name][index_name]=value
        
        return the_dict

    def generate_aitchison_pcoa(self, df):
        """
        Takes in count dataframe and generate pcoa matrix using aitchison distance matrix as input (aitchison distance is calculated in this function)
        
        Parameters
        ------------
        dataframe: pandas dataframe object
            where 
                rows = samples
                columns = OTU
                and each datapoint is the either read count or relative abundance. 

        Returns
        ------------
        pcoa.sample = pandas dataframe object
            this stores the pcoa generated values to be visualised
        
        dist_matrix = pandas dataframe object
            this stores the distance scores used to generate the pcoa values. 

        """
        dist_matrix = self.aitchison_distance_matrix(df)
        dm = DistanceMatrix(dist_matrix, dist_matrix.index)
        pcoa = ordination.pcoa(dm)
        
        return pcoa.samples, dist_matrix

    def get_kingdom_sets(self, dataframe):

        df = dataframe.copy()
        kingdomsets = {'virus': set(),
                'bacteria': set(),
                'eukaryote': set(),
                'archaea': set()
                }
        for taxid in df.columns:
            taxid = int(taxid)
            if 2 in ncbi.get_lineage(taxid):
                kingdomsets['bacteria'].add(taxid)
            elif 10239 in ncbi.get_lineage(taxid):
                kingdomsets['virus'].add(taxid)
            elif 2759 in ncbi.get_lineage(taxid):
                kingdomsets['eukaryote'].add(taxid)
            else:
                kingdomsets['archaea'].add(taxid)
            
        return kingdomsets


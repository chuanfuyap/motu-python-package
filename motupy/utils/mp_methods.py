"""
    This module contains methods for processing single sample read dictionary files.
"""
from ete3 import NCBITaxa
import pandas as pd
import numpy as np

ncbi = NCBITaxa
class utils():
    def taxa_cumulation(self, reads_dictionary, basic_ranks=['superkingdom',
                            'phylum',
                            'class',
                            'order',
                            'family',
                            'genus',
                            'species']):
        """
        Goes through all the taxa_id (keys) in the dict.
        For each of the taxa_id, it iteratively goes through their lineage using ete3 NCBITaxa.get_lineage().
        As it goes through lineage of taxa_id, it finds if there is a key in the dict corresponding to the taxa_id.
        If the key for the taxa_id is not found, it will be created and read counts from lower rank will be added to it.

        Parameters
        ------------
        reads_dictionary: dict,
            reads dictionary file such as the OTUdata

        basic_ranks: list [str],
            the desired ranks to be processed in the reads dictionary

        Returns
        ------------
        reads_dictionary: dict,
            reads dictionary file within class but now with taxa cumulation done. 

        """

        taxid_key_list = list(reads_dictionary.keys())
        for taxid_key in taxid_key_list:
            ordered_lineage = self.transform_lineage(taxid_key, basic_ranks)

            if taxid_key != ordered_lineage[0]: ## if the current taxid(taxid_key) is not the lowest rank (not strain/subsps i.e. anything beyond species level)
                reads_dictionary = self.update_dictfile(ordered_lineage[0], taxid_key, reads_dictionary, basic_ranks)

            for taxid in ordered_lineage[1:]:
                if taxid not in reads_dictionary:
                    reads_dictionary[taxid] = 0

        flipped_ranks = basic_ranks[::-1]
        ranks = self.update_ranks(reads_dictionary)

        for i, r in enumerate(flipped_ranks[:-1]): ## this reverses it to start from sps rank, so sps read count will be added to genus and genus to family and so on. 
            parent_rank = flipped_ranks[i+1]
            for taxid in ranks[r]:
                parent_taxid = self.find_parent_taxid(taxid, parent_rank, flipped_ranks)
                reads_dictionary[parent_taxid] += reads_dictionary[taxid]

        return reads_dictionary

    def transform_lineage(self, taxanomy_id, basic_ranks):
        """
        Takes a tax id and return a list of lineage of tax id ordered from lowest to highest.

        Parameters
        ------------
        taxanomy_id: int,
            NCBI Taxanomy ID

        basic_ranks: list [str],
            the desired ranks to be seen in lineage

        Returns
        ------------
        lineage: list,
            a list of linegae taxa_id ordered from lowest to highest taxanomy rank

        """
        
        tmp_lineage = []
        for taxid in ncbi.get_lineage(taxanomy_id):
            if ncbi.get_rank([taxid])[taxid] in basic_ranks:
                tmp_lineage.append(taxid) 

        return tmp_lineage[::-1]# invert list, so lowest rank appears first

    def update_ranks(self, reads_dictionary):
        """
        Process done right after reading otu data into a dictionary.
        This method collects all the taxa id in the dictionary's keys and group them into their respective taxanomy rank 
        and store them in a dictionary (self.ranks). 

        Parameters
        ------------
        reads_dictionary: dict,
            reads dictionary file such as the OTUdata

        Returns
        ------------
        ranks_group_dict: dict,
            the dictionary file where
                key = taxonomic rank
                value = list of taxa id belonging to the rank
        """
        ranks_group_dict = {'superkingdom': set(),
                            'phylum': set(),
                            'class': set(),
                            'order': set(),
                            'family': set(),
                            'genus': set(),
                            'species': set()}
        for taxid in reads_dictionary.keys():
            taxid = int(taxid)
            taxid_rank = ncbi.get_rank([taxid])[taxid] 
            if taxid_rank in ranks_group_dict:
                ranks_group_dict[taxid_rank].add(taxid)
        
        return ranks_group_dict

    def update_dictfile(self, taxanomy_id, dictionary_key, reads_dictionary, basic_ranks):
        """
        Takes a tax id and and update the dictinoary data by adding to existing key or make new one
        Parameters
        ------------
        taxanomy_id: int,
            NCBI Taxanomy ID from lineage list

        dictionary_key: int,
            NCBI Taxanomy ID from self.otufile dictionary

        reads_dictionary: dict,
            reads dictionary file such as the OTUdata

        basic_ranks: list [str],
            the desired ranks to be processed in the reads dictionary
        
        Returns
        ------------
        reads_dictionary: dict,
            reads dictionary file 

        """
        if taxanomy_id in reads_dictionary:
            reads_dictionary[taxanomy_id] += reads_dictionary[dictionary_key]

        else:
            reads_dictionary[taxanomy_id] = reads_dictionary[dictionary_key]
        
        taxid_rank = ncbi.get_rank([dictionary_key])[dictionary_key] 
        if taxid_rank not in basic_ranks:
            del reads_dictionary[dictionary_key]
        
        return reads_dictionary

    def find_parent_taxid(self, taxid, parent_rank, taxanomy_ranks):
        """
        Process done during taxa_cumulation step to identify the taxa if of the parent rank 
        E.g. given a taxa id of a species, it will find the taxa if of its genus. 

        Parameters
        ------------
        taxid: int, 
            the ncbi taxanomy id of the current taxa id being processed
        parent_rank: str,
            the parent rank in which to look for its taxa id
        taxanomy_ranks: list,
            the list of the 7 basic taxanomy ranks, this is used in the event the given taxa id does not have an immediate parent rank
            E.g. if a species does not have genus, it will look for family, if no family is found it will look for one rank above family, etc. until a rank is found. 

        Returns
        ------------
        parent_taxid: int,
            the ncbi taxa id of the rank above the input taxid

        """
        taxid_lineage = ncbi.get_lineage(taxid)
        rankID_dict = {v:k for k,v in ncbi.get_rank(taxid_lineage).items()} ## flipping taxid to rank dictionary to rank to taxid

        if parent_rank not in rankID_dict:
            index = taxanomy_ranks.index(parent_rank)
            index+=1    

            find_higher_rank = False
            while not find_higher_rank:
                higher_rank = taxanomy_ranks[index]
                if higher_rank in rankID_dict:
                    find_higher_rank = True
                    parent_rank = higher_rank
                index+=1
        try:
            parent_taxid = rankID_dict[parent_rank]
        except:
            print('[ERROR] find_parent_taxid:%s\t%s'%(parent_rank, rankID_dict))
        return parent_taxid

    def turn_reads_to_clade_relative_abundance(self, reads_dictionary, basic_ranks=['superkingdom',
                            'phylum',
                            'class',
                            'order',
                            'family',
                            'genus',
                            'species']):
        """
        Turns the read counts into clade relative abundance, i.e. a tax id a given rank will be divided by the total number of reads in the given rank
        E.g. a taxid in species rank with a given read count will be divided by total number of reads that makes up species rank. 

        Parameters
        ------------
        reads_dictionary: dict,
            reads dictionary file such as the OTUdata

        basic_ranks: list [str],
            the desired ranks to be processed in the reads dictionary

        Returns
        ------------
        N/A
        """
        ranks = self.update_ranks(reads_dictionary)
        clade_read_counts = self.get_clade_read_counts(reads_dictionary, basic_ranks, ranks)

        for rank in basic_ranks: ## to iterate through every single basic taxanomic ranks
            total_rank_reads = clade_read_counts[rank]
            for taxid in ranks[rank]: ## to iterate through every single taxid that belongs to the given taxanomic rank
                reads_dictionary[taxid] = reads_dictionary[taxid]/total_rank_reads

    def get_clade_read_counts(self, reads_dictionary, basic_ranks, ranks):
        """
        Generate a dictionary with total clade read counts, refer below:

        Parameters
        ------------
        reads_dictionary: dict,
            reads dictionary file such as the OTUdata

        basic_ranks: list [str],
            the desired ranks to be processed in the reads dictionary

        ranks: dict,
            ranks dictionary file where {key: a taxonomic rank, value: a set containing all the taxonomic ID belonging to the rank}

        Returns
        ------------
        clade_read_counts : dict,
            a dictionary where key = a taxanomy rank
                                value = total read counts in the given taxanomy rank. 
        """
        clade_read_counts = {}
        for rank in basic_ranks:
            for taxid in ranks[rank]:
                if rank in clade_read_counts:
                    clade_read_counts[rank] += reads_dictionary[taxid]
                else:
                    clade_read_counts[rank] = reads_dictionary[taxid]

        return clade_read_counts

    def select_rank(self, reads_dictionary, desired_ranks):
        """

        Parameters
        ------------
        reads_dictionary: dict,
            reads dictionary file such as the OTUdata

        desired_ranks: list [str],
            the desired ranks to be processed in the reads dictionary

        Returns
        ------------
        desired_rank_otufile : dict,
            the dict file containing taxa id keys from only the desired ranks
        """
        ranks = self.update_ranks(reads_dictionary)
        desired_rank_otufile = {}
        for rank in desired_ranks:
            rank = rank.lower()
            for taxid in ranks[rank]:
                desired_rank_otufile[taxid] = reads_dictionary[taxid]

        return desired_rank_otufile
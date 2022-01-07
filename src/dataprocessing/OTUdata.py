"""
    OTUdata object is for storing and processing single sample of OTU data.
    The data will be stored as a dictionary {key: OTU taxa id number, value: read counts/relative abundance}
"""
from .kaiju_output import *
from ete3 import NCBITaxa
from .replacement_taxa import get_dict
class OTUdata:
    def __init__(self, file_loc, input_type, artifact_threshold=0, verbose=0, extension=None, include_strains=False):
        """
        Creates OTUdata object for manipulation/transformation

        Parameters
        ------------
        file_loc: str,
            location/file name of the input file (taxa profiling output file)

        input_type: str,
            what type of taxa profiling tool was used to generate the file, "E.g. Kaiju"

        Returns
        ------------
        N/A

        """ 

        self.file_id = ''
        self.include_strains = include_strains
        self.verbose = verbose
        self.ncbi = NCBITaxa()
        self.cumulated = False
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

        if extension is None:
            self.file_id = self.process_file_name(file_loc)
        else:
            self.file_id = self.process_file_name_with_known_extension(file_loc, extension)

        if input_type.lower()  not in ['kaiju','old kaiju','old kaiju v2']:
            print("[ERROR] only 'kaiju','old kaiju','old kaiju v2' input type supported")
        if input_type.lower() == 'kaiju':
            self.otufile = read_kj(file_loc, artifact_threshold)
        
        if input_type.lower() == 'old kaiju':
            self.otufile = read_old_kj(file_loc, artifact_threshold)

        if input_type.lower() == 'old kaiju v2':
            self.otufile = read_old_kjV2(file_loc, artifact_threshold)

        #root_id = [-1, 1, 131567]
        root_id = [-1]
        for root in root_id:
            if root in self.otufile:
                del self.otufile[root]

        ### these two are problematic ncbi taxa id that have been deleted, I am now replacing them with their updated taxa id
        
        fix_taxa = get_dict()
        for k,v in fix_taxa.items():
            if k in self.otufile:
                self.otufile[v] = self.otufile[k]
                del self.otufile[k]

        self.update_ranks()
        self.update_superkingdom()

    
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
            the id_name to be associated with OTUdata object
        """
        slash_split = file_loc.split('/')
        slash_split = slash_split[-1]
        id_name = slash_split.replace(ext, '')
        
        return id_name


    def process_file_name(self, file_loc):
        """
        Parameters
        ------------
        file_loc : str,
            the name/location of the file being processed

        Returns
        ------------
        id_name: str,
            the id_name to be associated with OTUdata object
        """
        slash_split = file_loc.split('/')
        slash_split = slash_split[-1]
        dot_split = slash_split.split('.')
        
        id_name = dot_split[0]
        
        return id_name

    def taxa_cumulation(self):
        """
        Goes through all the taxa_id (keys) in the dict.
        For each of the taxa_id, it iteratively goes through their lineage using ete3 NCBITaxa.get_lineage().
        As it goes through lineage of taxa_id, it finds if there is a key in the dict corresponding to the taxa_id.
        If the key for the taxa_id is not found, it will be created and read counts from lower rank will be added to it.

        Parameters
        ------------
        self.otufile: dict,
            reads dictionary file within class 

        Returns
        ------------
        self.otufile: dict,
            reads dictionary file within class but now with taxa cumulation done. 

        """
        root_id = [-1, 1, 131567]

        for root in root_id:
            if root in self.otufile:
                del self.otufile[root]
        if not self.cumulated:
            taxid_key_list = list(self.otufile.keys())
            for taxid_key in taxid_key_list:
                ordered_lineage = self.transform_lineage(taxid_key)

                if taxid_key != ordered_lineage[0]: ## if the current taxid(taxid_key) is not the lowest rank (not strain/subsps i.e. anything beyond species level)
                    self.update_dictfile(ordered_lineage[0], taxid_key)

                    self.update_one_rank(ordered_lineage[0])
                    self.update_one_superkingdom(ordered_lineage[0])

                for taxid in ordered_lineage[1:]:
                    if taxid not in self.otufile:
                        self.otufile[taxid] = 0
                    
                    self.update_one_rank(taxid)
                    self.update_one_superkingdom(taxid)

            
            flipped_ranks = self.basic_ranks[::-1]

            for i, r in enumerate(flipped_ranks[:-1]): ## this reverses it to start from sps rank, so sps read count will be added to genus and genus to family and so on. 
                parent_rank = flipped_ranks[i+1]
                for taxid in self.ranks[r]:
                    parent_taxid = self.find_parent_taxid(taxid, parent_rank, flipped_ranks)
                    self.otufile[parent_taxid] += self.otufile[taxid]
                    
                    if self.verbose >=1 : print('[UPDATING PARENT ID READS] TAX ID %s \tNAME %s updated to %s with\nTAX ID %s \tNAME %s READS %s'%(parent_taxid, self.ncbi.get_taxid_translator([parent_taxid])[parent_taxid], self.otufile[parent_taxid], taxid, self.ncbi.get_taxid_translator([taxid])[taxid], self.otufile[taxid]))

            self.cumulated = True

    def transform_lineage(self, taxanomy_id):
        """
        Takes a tax id and return a list of lineage of tax id ordered from lowest to highest.
        Parameters
        ------------
        taxanomy_id: int,
            NCBI Taxanomy ID

        Returns
        ------------
        lineage: list,
            a list of linegae taxa_id ordered from lowest to highest taxanomy rank

        """
        tmp_lineage = []
        for taxid in self.ncbi.get_lineage(taxanomy_id):
            if self.ncbi.get_rank([taxid])[taxid] in self.basic_ranks:
                tmp_lineage.append(taxid) 

        return tmp_lineage[::-1]# invert list, so lowest rank appears first

    def update_dictfile(self, taxanomy_id, dictionary_key):
        """
        Takes a tax id and and update the dictinoary data by adding to existing key or make new one
        Parameters
        ------------
        taxanomy_id: int,
            NCBI Taxanomy ID from lineage list

        dictionary_key: int,
            NCBI Taxanomy ID from self.otufile dictionary
        
        Returns
        ------------
        N/A

        """
        if taxanomy_id in self.otufile:
            if self.verbose >= 2 : print('[ADDING READ COUNTS]TAX ID %s \tNAME %s already in dictionary with %s'%(taxanomy_id, self.ncbi.get_taxid_translator([taxanomy_id])[taxanomy_id], self.otufile[taxanomy_id]))

            self.otufile[taxanomy_id] += self.otufile[dictionary_key]

            if self.verbose >= 2 : print('TAX ID %s \tNAME %s updated to %s with\nTAX ID %s \tNAME %s READS %s'%(taxanomy_id, self.ncbi.get_taxid_translator([taxanomy_id])[taxanomy_id], self.otufile[taxanomy_id], dictionary_key, self.ncbi.get_taxid_translator([dictionary_key])[dictionary_key], self.otufile[dictionary_key]))

        else:
            if self.verbose >= 2 : print('[NEW ID ADDED] TAX ID %s \tNAME %s not found, adding to dictionary with\nTAX ID %s \tNAME %s READS %s'%(taxanomy_id, self.ncbi.get_taxid_translator([taxanomy_id])[taxanomy_id], dictionary_key, self.ncbi.get_taxid_translator([dictionary_key])[dictionary_key], self.otufile[dictionary_key]))

            self.otufile[taxanomy_id] = self.otufile[dictionary_key]
        
        taxid_rank = self.ncbi.get_rank([dictionary_key])[dictionary_key] 
        if taxid_rank not in self.basic_ranks:
            del self.otufile[dictionary_key]

            if self.verbose >= 1 : print('[NO RANK DELETION] TAX ID %s \tNAME %s deleted as it is not part of the desired ranks %s'%(dictionary_key, self.ncbi.get_taxid_translator([dictionary_key])[dictionary_key],taxid_rank.upper()))

    def update_superkingdom(self):
        """
        Process done right after reading otu data into a dictionary.
        This method collects all the taxa id in the dictionary's keys and group them into their respective superkingdom
        and store them in a dictionary (self.superkingdom). 

        Parameters
        ------------
        N/A

        Returns
        ------------
        N/A
        """
        for taxid in  self.otufile.keys():
            taxid = int(taxid)
            if 10239 in self.ncbi.get_lineage(taxid):
                self.superkingdom['virus'].add((taxid))
            elif 2 in self.ncbi.get_lineage(taxid):
                self.superkingdom['bacteria'].add((taxid))
            elif 2759 in self.ncbi.get_lineage(taxid):
                self.superkingdom['eukaryote'].add((taxid))
            elif 2157 in self.ncbi.get_lineage(taxid):
                self.superkingdom['archaea'].add((taxid))

    def update_one_superkingdom(self, taxid):
        """
        Process done during taxa_cumulation step as 'missing' parent OTUs (taxa id) are updated into the self.ranks dictionary

        Parameters
        ------------
        taxid: int,
            the ncbi taxanomy id to be added into self.ranks dictionary

        Returns
        ------------
        N/A
        """
        taxid = int(taxid)
        if 10239 in self.ncbi.get_lineage(taxid):
            self.superkingdom['virus'].add((taxid))
        elif 2 in self.ncbi.get_lineage(taxid):
            self.superkingdom['bacteria'].add((taxid))
        elif 2759 in self.ncbi.get_lineage(taxid):
            self.superkingdom['eukaryote'].add((taxid))
        elif 2157 in self.ncbi.get_lineage(taxid):
            self.superkingdom['archaea'].add((taxid))

    def update_ranks(self):
        """
        Process done right after reading otu data into a dictionary.
        This method collects all the taxa id in the dictionary's keys and group them into their respective taxanomy rank 
        and store them in a dictionary (self.ranks). 

        Parameters
        ------------
        N/A

        Returns
        ------------
        N/A
        """
        for taxid in self.otufile.keys():
            taxid = int(taxid)
            taxid_rank = self.ncbi.get_rank([taxid])[taxid] 
            if taxid_rank in self.ranks:
                self.ranks[taxid_rank].add(taxid)

    def update_one_rank(self, taxid):
        """
        Process done during taxa_cumulation step as 'missing' parent OTUs (taxa id) are updated into the self.ranks dictionary

        Parameters
        ------------
        taxid: int,
            the ncbi taxanomy id to be added into self.ranks dictionary

        Returns
        ------------
        N/A
        """
        taxid_rank = self.ncbi.get_rank([taxid])[taxid] 
        if taxid_rank in self.ranks:
            self.ranks[taxid_rank].add(taxid)

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
        taxid_lineage = self.ncbi.get_lineage(taxid)
        rankID_dict = {v:k for k,v in self.ncbi.get_rank(taxid_lineage).items()} ## flipping taxid to rank dictionary to rank to taxid

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

    def turn_reads_to_clade_relative_abundance(self):
        """
        Turns the read counts into clade relative abundance, i.e. a tax id a given rank will be divided by the total number of reads in the given rank
        E.g. a taxid in species rank with a given read count will be divided by total number of reads that makes up species rank. 

        Parameters
        ------------
        N/A

        Returns
        ------------
        N/A
        """

        if not self.cumulated:
            self.taxa_cumulation()

        clade_read_counts = self.get_clade_read_counts()

        for rank in self.basic_ranks: ## to iterate through every single basic taxanomic ranks
            total_rank_reads = clade_read_counts[rank]
            for taxid in self.ranks[rank]: ## to iterate through every single taxid that belongs to the given taxanomic rank
                self.otufile[taxid] = self.otufile[taxid]/total_rank_reads

    def get_clade_read_counts(self):
        """
        Generate a dictionary with total clade read counts, refer below:

        Parameters
        ------------
        N/A

        Returns
        ------------
        clade_read_counts : dict,
            a dictionary where key = a taxanomy rank
                                value = total read counts in the given taxanomy rank. 
        """
        clade_read_counts = {}
        for rank in self.basic_ranks:
            for taxid in self.ranks[rank]:
                if rank in clade_read_counts:
                    clade_read_counts[rank] += self.otufile[taxid]
                else:
                    clade_read_counts[rank] = self.otufile[taxid]

        return clade_read_counts

    def select_rank(self, desired_ranks):
        """

        Parameters
        ------------
        desired_ranks : list [str], 
            the desired taxanomy_ranks made up of str

        Returns
        ------------
        desired_rank_otufile : dict,
            the dict file containing taxa id keys from only the desired ranks
        """
        desired_rank_otufile = {}
        for rank in desired_ranks:
            rank = rank.lower()
            for taxid in self.ranks[rank]:
                desired_rank_otufile[taxid] = self.otufile[taxid]

        return desired_rank_otufile

    def identify_strain(self):
        """
        TO DO
        Parameters
        ------------
        N/A

        Returns
        ------------
        N/A
        """
        strain = False

        return strain
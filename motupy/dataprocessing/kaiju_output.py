"""
    Methods for importing kaiju summary output files
"""
def read_kj(filename, artifact_threshold=0):
    """
    Takes kaiju output file and make them into python dictionary.

    Parameters
    ------------
    filename: str,
        location/file name of the kaiju output file
    
    artifact_threshold: int,
        threshold for artifact range, if it is lower than threshold, the OTU is not added to the dictionary. 

    Returns
    ------------
    readsDict: dict,
        dictionary file where:
                key : NCBI taxanomy ID
                value : number of reads
    """ 
    readFile = open(filename, 'r')
    readsDict = {}
    all_lines = readFile.readlines()

    for line in all_lines[2:-3]:
        tokens = line.rstrip().split('\t')
        
        if len(tokens)==5:
            count = int(tokens[0])
            taxa_id = int(tokens[3])
            if count > artifact_threshold:
                readsDict[taxa_id]=int(count)

    return readsDict

def read_old_kj(filename, artifact_threshold=0):
    """
    Takes kaiju output file and make them into python dictionary.

    Parameters
    ------------
    filename: str,
        location/file name of the kaiju output file

    artifact_threshold: int,
        threshold for artifact range, if it is lower than threshold, the OTU is not added to the dictionary. 

    Returns
    ------------
    readsDict: dict,
        dictionary file where:
                key : NCBI taxanomy ID
                value : number of reads
    """ 
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    readFile = open(filename, 'r')
    readsDict = {}
    all_lines = readFile.readlines()

    for line in all_lines:
        tokens = line.rstrip().split('\t')
        
        if len(tokens)>2:
            count = int(tokens[0])
            final_id = -1 
            otu_name = tokens[final_id]
            while not bool(ncbi.get_name_translator([otu_name])):
                final_id-=1
                otu_name = tokens[final_id]
            
            taxa_id = ncbi.get_name_translator([otu_name])[otu_name]
            if count > artifact_threshold:
                readsDict[taxa_id[0]]=int(count)

    return readsDict

def read_old_kjV2(filename, artifact_threshold=0):
    """
    Takes kaiju output file and make them into python dictionary.

    Parameters
    ------------
    filename: str,
        location/file name of the kaiju output file

    artifact_threshold: int,
        threshold for artifact range, if it is lower than threshold, the OTU is not added to the dictionary. 

    Returns
    ------------
    readsDict: dict,
        dictionary file where:
                key : NCBI taxanomy ID
                value : number of reads
    """ 
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    readFile = open(filename, 'r')
    readsDict = {}
    all_lines = readFile.readlines()

    for line in all_lines[1:-3]:
        tokens = line.rstrip().split('\t')
        
        if len(tokens)==5:
            count = int(tokens[2])
            final_id = -1 

            taxa = tokens[4]
            taxa = taxa[:-1]
            taxa = taxa.split(';')
            taxa = [t.strip() for t in taxa]

            otu_name = taxa[final_id]
            while not bool(ncbi.get_name_translator([otu_name])):
                final_id-=1
                otu_name = taxa[final_id]
            
            taxa_id = ncbi.get_name_translator([otu_name])[otu_name]
            if count > artifact_threshold:
                readsDict[taxa_id[0]]=int(count)

    return readsDict
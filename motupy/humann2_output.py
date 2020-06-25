"""
    Methods for importing humann2 merged gene families output files
"""
def read_h2(filename):
    """
    Takes humann2 output file and make them into python dictionary.

    Parameters
    ------------
    filename: str,
        location/file name of the kaiju output file
    
    
    Returns
    ------------
    readsDict: dict,
        dictionary file where:
                key : uniprotID
                value : read counts?
    """ 
    readFile = open(filename, 'r')
    readsDict = {}
    all_lines = readFile.readlines()

    for line in all_lines[1:]:
        tokens = line.rstrip().split('\t')

        if "|" not in tokens[0]:
            readsDict[tokens[0]]=tokens[1]

    return readsDict
def bacteria_counter(reads_dictionary, ncbi):
    """
    counts the number of bacteria in the reads dictionary file

    Parameters
    ------------
    reads_dictionary: dict,
        reads dictionary file with key = taxa id, value = number of reads

    ncbi: NCBITaxa(),
        ncbi taxa tool from ete3

    Returns
    ------------
    bact_count: int,
        the number of reads that contributes to bacteria kingdom
    """
    bact_count = 0
    for key in reads_dictionary.keys():
        if 2 in ncbi.get_lineage(key):
            bact_count += reads_dictionary[key]
    return bact_count

def virus_counter(reads_dictionary, ncbi):
    """
    counts the number of viruses in the reads dictionary file

    Parameters
    ------------
    reads_dictionary: dict,
        reads dictionary file with key = taxa id, value = number of reads

    ncbi: NCBITaxa(),
        ncbi taxa tool from ete3

    Returns
    ------------
    virus_count: int,
        the number of reads that contributes to virus kingdom
    """
    virus_count = 0
    for key in reads_dictionary.keys():
        if 10239 in ncbi.get_lineage(key):
            virus_count += reads_dictionary[key]
    return virus_count

def eukaryota_counter(reads_dictionary, ncbi):
    """
    counts the number of eukaryotes in the reads dictionary file

    Parameters
    ------------
    reads_dictionary: dict,
        reads dictionary file with key = taxa id, value = number of reads

    ncbi: NCBITaxa(),
        ncbi taxa tool from ete3

    Returns
    ------------
    eukaryote_count: int,
        the number of reads that contributes to eukaryote kingdom
    """
    eukaryote_count = 0
    for key in reads_dictionary.keys():
        if 2759 in ncbi.get_lineage(key):
            eukaryote_count += reads_dictionary[key]
    return eukaryote_count

def archaea_counter(reads_dictionary, ncbi):
    """
    counts the number of archaea in the reads dictionary file

    Parameters
    ------------
    reads_dictionary: dict,
        reads dictionary file with key = taxa id, value = number of reads

    ncbi: NCBITaxa(),
        ncbi taxa tool from ete3

    Returns
    ------------
    archaea_count: int,
        the number of reads that contributes to archaea kingdom
    """
    archaea_count = 0
    for key in reads_dictionary.keys():
        if 2157 in ncbi.get_lineage(key):
            archaea_count += reads_dictionary[key]
    return archaea_count

def otu_counts(reads_dictionary, ncbi):
    """
    counts the number of OTUs that makes up the 4 main kingdoms. 

    Parameters
    ------------
    reads_dictionary: dict,
        reads dictionary file with key = taxa id, value = number of reads

    ncbi: NCBITaxa(),
        ncbi taxa tool from ete3

    Returns
    ------------
    vir_otu, bact_otu, euk_otu, archa_otu: int,
        the number of OTUs for the corresponding kingdoms
    """
    vir_otu = 0
    bact_otu = 0
    euk_otu = 0
    archa_otu = 0
    for key in reads_dictionary.keys():
        if 2 in ncbi.get_lineage(key):
            bact_otu+=1
        elif 10239 in ncbi.get_lineage(key):
            vir_otu+=1
        elif 2759 in ncbi.get_lineage(key):
            euk_otu+=1
        else:
            archa_otu+=1
    
    return vir_otu, bact_otu, euk_otu, archa_otu
    
def describe(reads_dictionary, cumulated, ncbi):
    """
    prints out the number of reads/OTUs and percentage of reads in total for the 4 superkingdoms:
    Viruses, Bacteria, Eukaryota, and Archaea

    Parameters
    ------------
    reads_dictionary: dict,
        reads dictionary file with key = taxa id, value = number of reads

    cumulated: boolean,
        if the reads_dictionary have been taxa_cumulated

    ncbi: NCBITaxa(),
        ncbi taxa tool from ete3

    Returns
    ------------
    N/A
    """
    if cumulated:
        try:
            virus_count = reads_dictionary[10239]
        except:
            virus_count = 0
        try:
            bacteria_count = reads_dictionary[2]
        except:
            bacteria_count = 0
        try:
            eukaryote_count = reads_dictionary[2759]
        except:
            eukaryote_count = 0
        try:
            archaea_count = reads_dictionary[2157]
        except:
            archaea_count = 0 

        total_reads = virus_count + bacteria_count + eukaryote_count + archaea_count
        v,b,e,a = otu_counts(reads_dictionary, ncbi)
        total_otu = v+b+e+a

        total_read_counts = 0
        for value in reads_dictionary.values():
            total_read_counts+=value

        print('Viruses:\t%s\t%.2f%% with %s of OTUs'%(virus_count,(virus_count/total_reads)*100, v))
        print('Bacteria:\t%s\t%.2f%% with %s of OTUs'%(bacteria_count,(bacteria_count/total_reads)*100, b))
        print('Eukaryota:\t%s\t%.2f%% with %s of OTUs'%(eukaryote_count,(eukaryote_count/total_reads)*100, e))
        print('Archaea:\t%s\t%.2f%% with %s of OTUs'%(archaea_count,(archaea_count/total_reads)*100, a))
        print('TOTAL Reads of the 4 main Kingdoms:\t%s' %total_reads )
        print('TOTAL OTUs:\t%s'%total_otu)
        print('Reads Cumulated:\t%s' %cumulated)
        print('TOTAL Reads after cumulation:\t%s'%total_read_counts)

    else:
        virus_count = virus_counter(reads_dictionary, ncbi)
        bacteria_count = bacteria_counter(reads_dictionary, ncbi)
        eukaryote_count = eukaryota_counter(reads_dictionary, ncbi)
        archaea_count = archaea_counter(reads_dictionary, ncbi)

        total_reads = virus_count + bacteria_count + eukaryote_count + archaea_count
        
        v,b,e,a = otu_counts(reads_dictionary, ncbi)
        total_otu = v+b+e+a
        print('Viruses:\t%s\t%.2f%% with %s of OTUs'%(virus_count,(virus_count/total_reads)*100, v))
        print('Bacteria:\t%s\t%.2f%% with %s of OTUs'%(bacteria_count,(bacteria_count/total_reads)*100, b))
        print('Eukaryota:\t%s\t%.2f%% with %s of OTUs'%(eukaryote_count,(eukaryote_count/total_reads)*100, e))
        print('Archaea:\t%s\t%.2f%% with %s of OTUs'%(archaea_count,(archaea_count/total_reads)*100, a))
        print('TOTAL Reads:\t%s' %total_reads)
        print('TOTAL OTUs:\t%s'%total_otu)
        print('Reads Cumulated:\t%s' %cumulated)


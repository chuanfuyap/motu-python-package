def makecombo(dnatraindf, dnatestdf, rnatraindf, rnatestdf, trainlabel, testlabel):
    ##make dna only data
    train, test = generateABT(dnatraindf, dnatestdf, trainlabel, testlabel)
    
    dna = {'train': train, 'test': test}
    ##make dna rna combined data
    dnatrain, dnatest = appendcolumn(dnatraindf, '_D'), appendcolumn(dnatestdf, '_D')
    rnatrain, rnatest = appendcolumn(rnatraindf, '_R'), appendcolumn(rnatestdf, '_R')
    
    dnatrain = dnatrain.reindex(rnatrain.index)
    dnatest = dnatest.reindex(rnatest.index)
    
    dnarnatrain = pd.concat([dnatrain, rnatrain])
    
    #dnarna = {'train': , 'test':}
    #active = {'train': , 'test':}

    return dna, dnatrain, rnatrain#, dnarna, active

def appendcolumn(dataframe, xna):
    df = dataframe.copy()
    df.columns = df.columns.map(lambda x: str(x)+xna)
    
    return df

def generateABT(traindf, testdf, trainlabel, testlabel):
    train = traindf.copy()
    test = testdf.copy()
    ##merging common OTUs
    merged = pd.concat([train, test], axis=0, sort=False)
    mergedfiltered = eda.filter_otu_coverage(merged, 0.1)
    
    train = mergedfiltered.loc[traindf.index]
    test = mergedfiltered.loc[testdf.index]
    
    train = matchCPlabels(train, trainlabel)
    
    test = pd.concat([test, testlabel], axis=1, sort=False)
    
    return train, test

def matchCPlabels(data, labels):
    df = data.copy()
    label = labels[['Subject', 'Asthma']].set_index('Subject')
    index = df.index
    
    match = {}
    
    for ind in index:
        i = int(ind[:3])
        match[ind] = label.loc[i].values[0]
        
    match = pd.Series(match, name='Asthma')
    
    df = pd.concat([df, match], axis=1, sort=False)

    return df
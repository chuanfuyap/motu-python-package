"""
    Timedecay analysis object for microbiome time series data.

    Method taken from this paper
    Shade, A., Gregory Caporaso, J., Handelsman, J. et al. A meta-analysis of changes in bacterial and archaeal communities with time. ISME J 7, 1493–1506 (2013). https://doi.org/10.1038/ismej.2013.54

    How: Log-linear model fitted to sample dissimilarities versus time between observations.
    Goal: Rate of community change (“turnover”)
    Pre-processing: Filtering of rare OTUs is recommended
    Remarks: Sensitive to time-scale

    NOTE: only dissimilarity measure implemented in this motupy package aitchison distance, hence the euclidean distance.
"""
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt 
import seaborn as sns; sns.set()

from skbio.stats import composition

from scipy.spatial.distance import pdist, euclidean
from itertools import combinations

import statsmodels.api as sm

from scipy import stats

class timedecay:
    def __init__(self, parameters=None):
        """
        Creates timedecay object to determine rate of community change over time
        i.e. if the microbiome gets more similar or more dissimilar over time. 

        Parameters
        ------------
        N/A
        
        Returns
        ------------
        N/A

        """
        self.parameters = parameters
        
    def fit(self, dataframe, metadataframe, verbose=True, inverse=False):
        """
        Performs timedecay analysis by fitting log-linear model
        
        Parameters
        ------------
        dataframe: pandas dataframe,
            microbiome count data
        metadataframe: pandas dataframe,
            metadata dataframe that should at least contain the subject column and time columns
            where subject column indicates which smaples belong to which subject 
                and the time column indicates the time of sampling. 
        
        Returns
        ------------
        lm: statsmodels LinearRegression fit results.
            the fitted linear regression model
        """

        abt = self.build_abt(dataframe, metadataframe,inverse)

        ## this adds an intercept to the linear model
        X = sm.add_constant(abt.timediff)
        y = abt.distance

        model = sm.OLS(y, X)
        results = model.fit()
        if verbose:
            print(results.summary())

        return results

    def fit_compare_intercept(self, dataframe, metadataframe, group="Group", labels=[1,0]):
        """
        DEPRECATED
        Fits data to for time decay analysis and compare intercept of two groups

        NOTE: function requires more testing to confirm that it is implemented correctly
        
        Parameters
        ------------
        dataframe: pandas dataframe,
            microbiome count data
        metadataframe: pandas dataframe,
            metadata dataframe that should at least contain the subject column and time columns
            where subject column indicates which smaples belong to which subject 
                and the time column indicates the time of sampling. 
            metadata must have 'Group' label
        group: str,
            the name of the column in metadataframe that contains group information
        labels: list
            made up of 2 elements containing the group labels, comparison is done in manner of labels[0] - labels[1]
        
        Returns
        ------------
        lm: statsmodels LinearRegression fit results.
            the fitted linear regression model
        """
        groupAix = metadataframe[metadataframe["Group"]==labels[0]].index
        groupBix = metadataframe[metadataframe["Group"]==labels[1]].index

        abtA = self.build_abt(dataframe.loc[groupAix], metadataframe.loc[groupAix])
        abtB = self.build_abt(dataframe.loc[groupBix], metadataframe.loc[groupBix])

        groupA = sm.add_constant(abtA.timediff)
        groupB = sm.add_constant(abtB.timediff)

        groupB["const"]=0

        X = pd.concat([groupA, groupB])
        X.columns = ["group", "timediff"]

        X = sm.add_constant(X)
        y = pd.concat([abtA.distance, abtB.distance])

        model = sm.OLS(y, X)
        results = model.fit()
        print(results.summary())

        return results

    def fit_compare_slope(self, dataframe, metadataframe, group="Group", labels=[1,0]):
        """
        DEPRECATED
        Fits data to for time decay analysis and compare intercept of two groups
        
        Parameters
        ------------
        dataframe: pandas dataframe,
            microbiome count data
        metadataframe: pandas dataframe,
            metadata dataframe that should at least contain the subject column and time columns
            where subject column indicates which smaples belong to which subject 
                and the time column indicates the time of sampling. 
            metadata must have 'Group' label
        group: str,
            the name of the column in metadataframe that contains group information
        labels: list
            made up of 2 elements containing the group labels, comparison is done in manner of labels[0] - labels[1]
        
        Returns
        ------------
        lm: statsmodels LinearRegression fit results.
            the fitted linear regression model
        """
        groupAix = metadataframe[metadataframe["Group"]==labels[0]].index
        groupBix = metadataframe[metadataframe["Group"]==labels[1]].index

        abtA = self.build_abt(dataframe.loc[groupAix], metadataframe.loc[groupAix])
        abtB = self.build_abt(dataframe.loc[groupBix], metadataframe.loc[groupBix])

        groupA = sm.add_constant(abtA.timediff)
        groupB = sm.add_constant(abtB.timediff)

        groupB["const"]=0

        X = pd.concat([groupA, groupB])
        X.columns = ["group", "timediff"]

        X = sm.add_constant(X)
        X["group*time"] = X.group*X.timediff
        y = pd.concat([abtA.distance, abtB.distance])

        model = sm.OLS(y, X)
        results = model.fit()
        print(results.summary())

        return results

    def build_abt(self, dataframe, metadataframe, inverse=False):
        """
        Builds the analytical base table needed for log linear model fitting
        
        Parameters
        ------------
        dataframe: pandas dataframe,
            microbiome count data
        metadataframe: pandas dataframe,
            metadata dataframe that should at least contain the subject column and time columns
            where subject column indicates which smaples belong to which subject 
                and the time column indicates the time of sampling. 
        
        Returns
        ------------
        xydf: pandas dtaframe,
            the analytical base table needed for time-decay model. 
        """

        df = dataframe.copy()
        df = self.clrtransform(df) ##clr transform count data
        
        metadata = metadataframe.copy()
        subjects = metadataframe.Subject.unique()
        
        xylist = []
        
        ## calculates sample distance within each subjects
        for i in subjects:
            ix = metadata[metadata.Subject == i].index
            subdf = df.loc[ix]
            submeta = metadata.loc[ix, "Time"]
            
            XY = self.scoreXY(subdf, submeta)
            XY["subject"] = [i]*XY.shape[0]
            
            xylist.append(XY)
        
        # joins the distance table from all subjects
        xydf = pd.concat(xylist).reset_index().drop(columns="index")
        
        if inverse:
            xydf["distance"] = 1/xydf["distance"]
        xydf["distance"] = np.log(xydf["distance"])
        
        # removes values of zero distance
        xydf=xydf[~np.isinf(xydf)]
        xydf.dropna(axis=0, how='any', inplace=True)

        return xydf

    def build_abt_distance(self, dataframe, metadataframe):
        """
        Builds the analytical base table needed for log linear model fitting
        
        Parameters
        ------------
        dataframe: pandas dataframe,
            microbiome count data
        metadataframe: pandas dataframe,
            metadata dataframe that should at least contain the subject column and time columns
            where subject column indicates which smaples belong to which subject 
                and the time column indicates the time of sampling. 
        
        Returns
        ------------
        xydf: pandas dtaframe,
            the analytical base table needed for time-decay model. 
        """

        df = dataframe.copy()
        df = self.clrtransform(df) ##clr transform count data
        
        metadata = metadataframe.copy()
        subjects = metadataframe.Subject.unique()
        
        xylist = []
        
        ## calculates sample distance within each subjects
        for i in subjects:
            ix = metadata[metadata.Subject == i].index
            subdf = df.loc[ix]
            submeta = metadata.loc[ix, "Time"]
            
            XY = self.scoreXY(subdf, submeta)
            XY["subject"] = [i]*XY.shape[0]
            
            xylist.append(XY)
        
        # joins the distance table from all subjects
        xydf = pd.concat(xylist).reset_index().drop(columns="index")
        
        xydf["distance"] = np.log(xydf["distance"])
        
        # removes values of zero distance
        xydf=xydf[~np.isinf(xydf)]
        xydf.dropna(axis=0, how='any', inplace=True)
        
        return xydf
    
    def clrtransform(self, dataframe):
        """
        Performs multiplicative replacement to fill in the zeroes followed by centred-log-ratio (clr) transformation.
        
        Parameters
        ------------
        dataframe: pandas dataframe,
            microbiome count data

        
        Returns
        ------------
        X_clr: pandas dataframe,
            dataframe containing the clr transformed values of the count data. 

        """
        df = dataframe.copy()    
        ### multiplicative replacement to make sure no zero
        df.fillna(0, inplace=True)
        X_imputed = composition.multiplicative_replacement(df)
        ### clr transform data
        X_clr = composition.clr(X_imputed)

        return pd.DataFrame(X_clr, columns=df.columns, index=df.index).sort_index()

    def scoreXY(self,subdataframe, subtime):
        """
        Generates the distance (Y) and time difference (X) measure needed for the log-linear model
        
        Parameters
        ------------
        subdataframe: pandas dataframe,
            subset of the count dataframe containing only the samples from a subject across time

        subtime: pandas dataframe,
            subset of the metadata datframe containing only the time of sampling for the subject's samples
        
        Returns
        ------------
        xy: pandas datframe,
            the dataframe with timediff between samples and distance between samples

        """
        ## X is the time difference (timediff)
        ## Y is the log value of sample distance (distance)
        df = subdataframe.T.copy()

        colnames = df.columns
        xy = pd.DataFrame(columns=[ "timediff", "distance"])

        for col1, col2 in combinations(colnames, 2):
            dist = euclidean(df[col1], df[col2])
            time = subtime[col2] - subtime[col1]
            xy = xy.append({"timediff":time, "distance":dist}, ignore_index=True)

        return xy

    def group_comparison_Ttest(self, dataset, metadata, group,inverse=False):
        """
        Performs T test to compute the statistical significance on the difference of time decay between two groups.
        This is done by looping through each subject within a group and computing their decay rate, and finally comparing decay rate of the two groups. 
        
        Parameters
        ------------
        dataset: pandas dataframe,
            the compositional data where,
                row = samples
                columns = OTUs

        metadata: pandas dataframe,
            the metadata of the file, where the row index must be shared with the dataset. 
            the columns will 
        
        group: str,
            the group within the metadata columns which will be compared 

        Returns
        ------------
        t: float, 
            t-score
        p: float,
            p-value
        decayrateA, decayrateB: list of float,
            the computed decay rate mean for group A and B
        interceptA, interceptB: list of float, 
            the computed intercept mean for group A and B

        """
        A, B = metadata[group].unique()
        groupAix = metadata[metadata[group]==A].index
        groupBix = metadata[metadata[group]==B].index
        
        abtA = self.build_abt(dataset.loc[groupAix], metadata.loc[groupAix],inverse)
        abtB = self.build_abt(dataset.loc[groupBix], metadata.loc[groupBix],inverse)
        
        decayrateA = []
        interceptA = []
        
        for sub in abtA.subject.unique():
            subABT = abtA[abtA.subject==sub]
            
            X = sm.add_constant(subABT.timediff)
            y = subABT.distance
            
            results = sm.OLS(y, X).fit()
            decayrate = results.params['timediff']
            decayrateA.append(decayrate)
            interceptA.append(results.params['const'])
            
        decayrateB = []
        interceptB = []
        
        for sub in abtB.subject.unique():
            subABT = abtB[abtB.subject==sub]
            
            X = sm.add_constant(subABT.timediff)
            y = subABT.distance
            
            results = sm.OLS(y, X).fit()
            decayrate = results.params['timediff']
            decayrateB.append(decayrate)
            interceptB.append(results.params['const'])
            
        t,p = stats.ttest_ind(decayrateA,decayrateB)
        
        return t, p, decayrateA, decayrateB, interceptA, interceptB
"""
    This module contains methods for Visualisation of analysis. 

    At the moment is it mostly for visualising aitchison distance comparison between samples. 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns; sns.set()

from scipy import stats

from itertools import combinations

class visualise():
    def within_between_group(self, distance_matrix, metadata, group):
        """
        compiles distances within the selected group and the distances between groups 
        
        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        group : str,
                    group/column within the metadata. 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the comparison within and between groups. 
        """
        within_groups = {}
        between_groups = []
        
        metalabels = metadata[group].unique()
        
        for label in metalabels:
            label_index = metadata[group]==label
            if distance_matrix.loc[label_index].shape[0]>0:
                lower_tri = np.tril(distance_matrix.loc[label_index,label_index])
                within_groups['Within %s' %label] = lower_tri[np.nonzero(lower_tri)]

        
        for combo in combinations(metalabels,2):
            labels1 = metadata[group]== combo[0]
            labels2 = metadata[group]== combo[1]
            
            between_df = distance_matrix.loc[labels1,labels2]
            between_groups.extend(between_df.values.flatten())
        
        within_groups['Between %s' %group] = between_groups
        
        tmp_df = []
        for k,v in within_groups.items():
            tmp_df.append(pd.DataFrame(v, columns=[k]))
            
        df = pd.concat(tmp_df, axis=1)
        
        return df

    def plot_violin_within_between_group(self, distance_matrix, metadata, group, title, savefile):
        """
        Generates the violin plot for within selected groups and between selected groups

        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        group : str,
                    group/column within the metadata. 

        title : str,
                    title of the plot generated
        
        savefile : str,
                    name/location of the plot to save as 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the comparison within and between groups. 
        """
        dataframe = self.within_between_group(distance_matrix, metadata, group)
        plt.figure(figsize=(12,7))
        sns.violinplot(data=dataframe)
        plt.title(title, fontsize=20)
        plt.ylabel('Aitchison Distance')
        plt.savefig(savefile)
        plt.show()

        return dataframe

    def within_between_group_mean(self, distance_matrix, metadata, group):
        """
        Compiles the within and between group distance comparison and compute a mean for each of them.
        
        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        group : str,
                    group/column within the metadata. 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        within_groups = []
        between_groups = []
        data = {}
        
        metalabels = metadata[group].unique()
        
        for label in metalabels:
            label_index = metadata[group]==label
            if distance_matrix.loc[label_index].shape[0]>0:
                lower_tri = np.tril(distance_matrix.loc[label_index,label_index])
                within_groups.append(lower_tri[np.nonzero(lower_tri)].mean())

        
        for combo in combinations(metalabels,2):
            labels1 = metadata[group]== combo[0]
            labels2 = metadata[group]== combo[1]
            
            between_df = distance_matrix.loc[labels1,labels2]
            between_groups.append(between_df.values.flatten().mean())
            
            
        data['Within'] = within_groups
        data['Between'] = between_groups
        
        tmp_df = []
        for k,v in data.items():
            tmp_df.append(pd.DataFrame(v, columns=[k]))
            
        df = pd.concat(tmp_df, axis=1)
        
        return df

    def barplot_within_between_categories(self, distance_matrix, metadata, categories, title, savefile):
        """
        Generates barplots for mean of within/between group comparison.

        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        categories : list,
                    list of str of all the categories to be plotted for the barplot

        title : str,
                    title of the plot generated
        
        savefile : str,
                    name/location of the plot to save as 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        df_list =[]
        for c in categories:
            df = self.within_between_group_mean(distance_matrix, metadata, c)
            melted = pd.melt(df)
            melted.columns = ['Legend', 'Aitchison Distance']
            melted['Category'] = [c]*melted.shape[0]
        
            df_list.append(melted)
            
        combined_df = pd.concat(df_list, axis=0)
        combined_df.dropna(inplace=True)
        
        plt.figure(figsize=(12,7))
        sns.barplot(data=combined_df, x='Category', y='Aitchison Distance', hue='Legend')
        plt.title(title, fontsize=20)
        plt.savefig(savefile)
        plt.show()
        
        return combined_df

    def withinSubjectsGroups(self, distance_matrix, metadata, group, label):
        """
        Compiles for mean/std of within subject comparisons within a group
        
        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        group : str,
                    group/column within the metadata. 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        subjects = metadata[metadata[group]==label].Subject.unique()
        
        subject_distance_array = {}
        for sub in subjects:
            sub_index = metadata.Subject==sub
            lower_tri = np.tril(distance_matrix.loc[sub_index,sub_index])
            subject_distance_array[sub] = lower_tri[np.nonzero(lower_tri)]

        subject_distance_df = pd.DataFrame(columns=subject_distance_array.keys())
        for col in subject_distance_df.columns:
            subject_distance_df[col]=pd.Series(subject_distance_array[col])
            
        mean_array = subject_distance_df.mean(0) 
        std_array = subject_distance_df.std(0)
            
        return mean_array.values, std_array.values#, subject_distance_df

    def violin_group_meanstd(self, distance_matrix, metadata, group, title, savefile):
        """
        Generates violin plots for mean/std of within subject comparisons within a group

        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        categories : list,
                    list of str of all the categories to be plotted for the barplot

        title : str,
                    title of the plot generated
        
        savefile : str,
                    name/location of the plot to save as 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        meandf = {}
        stddf = {}
        
        for label in metadata[group].unique():
            meanseries, stdseries = self.withinSubjectsGroups(distance_matrix, metadata, group, label)
            
            meandf[label] = meanseries
            stddf[label] = stdseries
            
        meandf = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in meandf.items() ]))
        stddf = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in stddf.items() ]))
        
        plt.figure(figsize=(12,7))
        sns.violinplot(data=meandf)
        plt.title("%s mean"% title, fontsize=20)
        plt.ylabel('Aitchison Distance')
        plt.savefig("%s_mean.png"% savefile)
        plt.show()
        
        plt.figure(figsize=(12,7))
        sns.violinplot(data=stddf)
        plt.title("%s standard deviation"% title, fontsize=20)
        plt.ylabel('Aitchison Distance')
        plt.savefig("%s_std.png"% savefile)
        plt.show()
        
        return meandf, stddf

    def within_between_subjects(self, distance_matrix, metadata, group):
        """
        Compiles for comparison each within subject and between subjects. 
        
        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        group : str,
                    group/column within the metadata. 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        within_groups = {}
        between_groups = []
        
        metalabels = metadata[group].unique()
        
        for label in metalabels:
            label_index = metadata[group]==label
            if distance_matrix.loc[label_index].shape[0]>0:
                lower_tri = np.tril(distance_matrix.loc[label_index,label_index])
                within_groups['Within %s' %label] = lower_tri[np.nonzero(lower_tri)]

        
        for combo in combinations(metalabels,2):
            labels1 = metadata[group]== combo[0]
            labels2 = metadata[group]== combo[1]
            
            between_df = distance_matrix.loc[labels1,labels2]
            between_groups.extend(between_df.values.flatten())
        
        within_groups['Between %s' %group] = between_groups
        
        tmp_df = []
        for k,v in within_groups.items():
            tmp_df.append(pd.DataFrame(v, columns=[k]))
            
        df = pd.concat(tmp_df, axis=1)
        mean_values = df.mean(0)
        mean_values.drop('Between %s' %group, inplace=True)
        sorted_columns = mean_values.sort_values().index.append(pd.Index(['Between Subject']))
        
        
        return df[sorted_columns]

    def plot_violin_within_subjects(self, distance_matrix, metadata, group, title, savefile):
        """
        Generates violin plots for within each subject and between subjects comparison. 

        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        categories : list,
                    list of str of all the categories to be plotted for the barplot

        title : str,
                    title of the plot generated
        
        savefile : str,
                    name/location of the plot to save as 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        dataframe = self.within_between_subjects(distance_matrix, metadata, group)
        plt.figure(figsize=(21,10))
        sns.violinplot(data=dataframe)
        plt.title(title, fontsize=20)
        plt.xticks(rotation=90)
        plt.ylabel('Aitchison Distance')
        plt.savefig(savefile)
        plt.show()
        
        return dataframe

    def within_between_groupV3(self, distance_matrix, metadata, group):
        """
        Compiles for mean of each within subject and between subjects. 
        
        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        group : str,
                    group/column within the metadata. 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        datadict = {}
        within_groups = []
        between_groups = []
        
        metalabels = metadata[group].unique()
        
        for label in metalabels:
            label_index = metadata[group]==label
            if distance_matrix.loc[label_index].shape[0]>0:
                lower_tri = np.tril(distance_matrix.loc[label_index,label_index])
                lower_tri = lower_tri[np.nonzero(lower_tri)]
                within_groups.append(lower_tri.flatten().mean())
                
        
        datadict['Within %s' %group] = within_groups

        
        for combo in combinations(metalabels,2):
            labels1 = metadata[group]== combo[0]
            labels2 = metadata[group]== combo[1]
            
            between_df = distance_matrix.loc[labels1,labels2]
            between_groups.append(between_df.values.flatten().mean())
        
        datadict['Between %s' %group] = between_groups
        
        tmp_df = []
        for k,v in datadict.items():
            tmp_df.append(pd.DataFrame(v, columns=[k]))
            
        df = pd.concat(tmp_df, axis=1)
        
        
        return df

    def plot_violin_within_betweenV3(self, distance_matrix, metadata, group, title, savefile):
        """
        Generates violin plots for mean of within each subject and between subjects comparison. 

        Parameters
        ------------
        distance_matrix : pandas DataFrame,
                        aithcison distance matrix generated from aitchison_distance_matrix method in EDA module. 

        metadata : pandas DataFrame,
                        the metadata for the sample highlighting the group of each samples
                        rows = samples
                        columns = grouping info, (e.g. Treatment Groups, labeled accoridingly)

        categories : list,
                    list of str of all the categories to be plotted for the barplot

        title : str,
                    title of the plot generated
        
        savefile : str,
                    name/location of the plot to save as 

        Returns
        ------------
        dataframe : pandas DataFrame generated following the mean calculation within and between groups. 
        """
        dataframe = self.within_between_groupV3(distance_matrix, metadata, group)
        plt.figure(figsize=(12,7))
        sns.violinplot(data=dataframe)
        plt.title(title, fontsize=20)
        plt.ylabel('Mean Aitchison Distance')
        plt.savefig(savefile)
        plt.show()
        
        return dataframe 


    def mean_confidence_interval(self, data, confidence=0.95):
        """
        Compute confidence interval given a list of values
        
        Parameters
        ------------
        data: list of values

        confidence: float,
            fractional value to compute confidence interval in

        Returns
        ------------
        m: float, 
            mean value
        m-h: float,
            lower confidence interval value
        m+h: float,
            upper confidence interval value

        """
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), stats.sem(a)
        h = se * stats.t.ppf((1 + confidence) / 2., n-1)
        return m, m-h, m+h

    def plotvalues(self, slope, intercept):
        """
        Compute a straight line plot given the slope and intercept values
        
        Parameters
        ------------
        slope: float,

        intercept: float,

        Returns
        ------------
        xaxis: list of float,
        yaxis: list of float

        """
        xaxis = np.arange(0,29)
        
        yaxis = intercept+slope*xaxis
        
        return xaxis, yaxis


    def plot_timedecay(self, slopeA, interceptA, slopeB, interceptB, dataA, dataB, labelA, labelB, pvalue):
        """
        Generates regression plot by computing mean value of a list of slope and intercept values for two groups
        
        Parameters
        ------------
        slopeA, slopeB: list of float,

        interceptA, interceptB: list of float,

        dataA, dataB: pandas DataFrame,
            the abt used in computing the decay rates, used for scatter plot accompanying the regression plot.
            generated using build_abt function. 

        labelA, labelB: str,
            labels for the Groups to be used in the plot.
        
        pvalue: float,
            pvalue computed using group_comparison_Ttest to be presented on the title of the plot.

        Returns
        ------------
        N/A

        """
        
        slope, slopedown, slopeup = self.mean_confidence_interval(slopeA)
        intercept, interceptdown, interceptup = self.mean_confidence_interval(interceptA)
        
        xA, yA = self.plotvalues(slope, intercept)
        _, yAdown = self.plotvalues(slopedown,interceptdown)
        _, yAup = self.plotvalues(slopeup,interceptup)

        slope, slopedown, slopeup = self.mean_confidence_interval(slopeB)
        intercept, interceptdown, interceptup = self.mean_confidence_interval(interceptB)
        
        xB, yB = self.plotvalues(slope, intercept)
        _, yBdown = self.plotvalues(slopedown,interceptdown)
        _, yBup = self.plotvalues(slopeup,interceptup)
        
        
        fig, (density, straightline) = plt.subplots(1, 2,figsize=(15,7),gridspec_kw={'width_ratios': [1, 2]})

        fig.suptitle("Time Decay, pvalue={}".format(pvalue),fontsize=15)

        straightline.plot(xA, yA, color='red', label=labelA)
        straightline.fill_between(xA, (yAdown), (yAup), alpha=.1, color='red')

        straightline.plot(xB, yB, label=labelB)        
        straightline.fill_between(xB, (yBdown), (yBup), color='b', alpha=.1)
        
        sns.scatterplot(x='timediff', y='distance', color='red',data=dataA, ax=straightline)
        sns.scatterplot(x='timediff', y='distance', color='blue',data=dataB, ax=straightline)
        
        straightline.set_xlabel("Time difference")
        straightline.set_ylabel("Aitchison Distance (log)")

        straightline.set_title("Time Decay plot")

        sns.distplot(slopeA, color='r', label=labelA, ax=density)
        sns.distplot(slopeB, color='b', label=labelB, ax=density)

        density.set_title("Decay Histogram/Distribution Plot")
        density.set_xlabel("Decay rate")

        plt.legend()
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
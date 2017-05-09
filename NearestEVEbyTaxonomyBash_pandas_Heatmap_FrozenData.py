# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 16:23:56 2016

@author: zwhitfield
Takes output(statistical enrichment for TEs overlapping EVEs derived from a specific viral family/genus/species etc...) from
 NearestEVEbyTaxonomyBash_pands_onlyOverlap_Functions_FrozenData.py
 Creates heatmap based on binomial test pValue enrichment using Seaborn (could also use Fishers exact test pValue. Numbers are always similar).
"""
import sys
import pandas as pd
import numpy as np
import scipy.stats as sp
#scipy.__version__
#np.__version__
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
#sns.__version__
#sns.set(context="paper", font="monospace")

#matplotlib.__version__
plt.style.use('ggplot')

inputdir = str(sys.argv[1])
outputdir = str(sys.argv[2])

#This will be used to only select specific entries in dataset with TEs of the given type.
#Can be at any taxonomy level (class, subclass, family, element, etc...)
#The specific level is specified by filteredBy_EVEs (ie name of column to subset using filteredBy)
filteredBy_TEs = str(sys.argv[3])
currentClassification = str(sys.argv[4])
filteredBy_EVEs = str(sys.argv[5])
analysisType = str(sys.argv[6])

# inputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/forStats/"
# outputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/heatmaps/"
# filteredBy_TEs="NONE" #Enter 'NONE' if don't want any filtering.
# currentClassification="TEclass"
# filteredBy_EVEs="family"
# analysisType = "overlapOrNearest" # "overlapOrNearest" OR "nearestOnly" OR "overlapOnly"

#possibleOrientations = ['same','opposite']
possibleOrientations = ['same']

for currentOrientation in possibleOrientations:
    allTEsOfInterest = []
    heatmapProportionCutoff = 0.1 #This is cutoff for displaying a category. If TE particular element is less than this proportion in all viral categories (virusFamilies), it is not shown

    #Specify what EVE-taxonomy labels to put on heatmap. Variable is called virusFamilies, but can be at the level of
    #species, family, group etc...
    # virusFamilies = ['Xincheng anphevirus', 'Quang Binh virus',
    #        'Grass carp rhabdovirus V76', 'Wuhan Mosquito Virus 8',
    #        'Wuhan Mosquito Virus 9', 'Yongjia Tick Virus 2',
    #        'Nkolbisson virus', 'Kamiti River virus', 'PFRV', 'Cilv-C',
    #        'Mosqueiro virus', 'Wutai Mosquito Virus', 'Le Dantec virus',
    #        'Tench rhabdovirus S64', 'Scophthalmus maximus rhabdovirus',
    #        'Nienokoue virus', 'CPV', 'Malakal virus', 'Berrimah virus', 'CJVS',
    #        'CFA', 'Kotonkan virus', 'Marco virus', 'OBOV',
    #        'Wuchang Cockraoch Virus 3', 'Aedes flavivirus',
    #        'Culex tritaeniorhynchus rhabdovirus', 'Tupaia rhabdovirus',
    #        'Barur virus', 'Blueberry necrotic ring blotch virus',
    #        'Tibrogargan virus', 'Garba virus',
    #        'Grapevine leafroll-associated virus 7', 'Wuhan Mosquito Virus 2',
    #        'Oita virus']
           
    #virusFamilies = ['Rhabdoviridae', 'Flaviviridae', 'Closteroviridae', 'Chuviridae', 'Bunyaviridae', 'Virgaviridae']
    virusFamilies = ['Rhabdoviridae', 'Flaviviridae', 'Chuviridae']
    
    #First get list of all TE elements associated with EVEs of given family (ie they have enrichment scores assigned to them)
    #Creates single-row dataframe, with columns being the viral categories of interest. Will hold the names of all
    #TE elements (categorized under currentClassification) overlapping EVEs derived from each viral category
    df_TEsPerEVE = pd.DataFrame(index = [0], columns=virusFamilies)    
    for virusFamily in virusFamilies:
        enrichmentScores = pd.read_csv(inputdir + 'ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy_TEs + 'AND' + virusFamily  + '_StatsComparedToWholeGenome_' + currentOrientation + 'Strand_' + analysisType + '.txt', sep="\t")
        df_TEsPerEVE.ix[0][virusFamily] = list(enrichmentScores[currentClassification])
        allTEsOfInterest.append(list(enrichmentScores[currentClassification]))
    
    # Create single list of unique TE elements near ALL EVE families of interest
    allTEsOfInterest = [item for sublist in allTEsOfInterest for item in sublist]    
    allTEsOfInterest = list(pd.unique(allTEsOfInterest))
    
    #creat initial dataframe of all possible TE elements (columns) and all virus families of interest (rows)
    #Fill with 0.0 as initial proportion, will be replaced later with proportion of given TE element compared to rest of TE elements found near EVEs derived from given viral family.
    df_values = pd.DataFrame(index=[virusFamilies], columns=allTEsOfInterest)
    df_values = df_values.fillna(float(0.0))
    
    #create initial dataframe of whether enrichment is significant.
    #The ' ' will be replaced with '*' if pvalue < than 0.0001
    df_significance = pd.DataFrame(index=[virusFamilies], columns=allTEsOfInterest)
    df_significance = df_significance.fillna(' ')
    
    #create initial dataframe of whether a given TE element was in the dataset for a given EVE viral family
    #will be changed to False if TE IS present. This value is eventually passed to seaborn about whether to mask it or not.
    df_present = pd.DataFrame(index=[virusFamilies], columns=allTEsOfInterest)
    df_present = df_present.fillna(True)
    
    
    #Fill in above dataframes with relevant/actual information.
    for virusFamily in virusFamilies:
        enrichmentScores = pd.read_csv(inputdir + 'ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy_TEs + 'AND' + virusFamily  + '_StatsComparedToWholeGenome_' + currentOrientation + 'Strand_' + analysisType + '.txt', sep="\t")
        #Loop through each listed TE element in 'enrichmentScores'
        for row in range(0,len(enrichmentScores[currentClassification])):
            currentElement = enrichmentScores[currentClassification][row]
            currentSignificanceScore = enrichmentScores['pValueOneSidedBinom'][row] #By binomial test p-value
            if currentElement in df_TEsPerEVE.ix[0][virusFamily]: #this is a redundant line I think, but good to check
                df_present.ix[virusFamily][currentElement] = False
                if currentSignificanceScore <= 0.0001: #pValue was arbitrarily chosen
                    df_significance.ix[virusFamily][currentElement] = '*'
    
            currentProportion = float(enrichmentScores['specificTEcounts'][row])/float(enrichmentScores['totalTEcountsOfType_'+ filteredBy_TEs][row]) #By proportion        
            df_values.ix[virusFamily][currentElement] = currentProportion

    # Remove TE elements (ie columns) from datasets if their proportion was <heatmapProportionCutoff for all virus families
    for TEelement in allTEsOfInterest:
        if len(df_values[(df_values[TEelement]<heatmapProportionCutoff)])==len(virusFamilies):
            df_values = df_values.drop([TEelement],1)
            df_significance = df_significance.drop([TEelement],1)
            df_present = df_present.drop([TEelement],1)
    
    #Generate heatmap, values based on proportion
    #A '*' is placed when binomial test p-value was significant
    #Will be grey if TE was not present in dataset for that viral family
    fig = plt.figure(figsize = (10, 10), facecolor='white')
    #with sns.axes_style("white"):
    sns.heatmap(df_values, 
                cmap = "YlGnBu", 
                #vmax=1,  
                annot = df_significance, 
                annot_kws={"size": 10},
                fmt = 's', 
                mask = df_present, #use this to grey out boxes of TE elements that don't appear at all for a given EVE viral classification
                square = True,
                linewidths=.5)
    #fig.set_size_inches(9,7)
    plt.yticks(rotation='horizontal')
    plt.xticks(rotation='vertical')
    fig.savefig(outputdir + 'ClassifiedBy_' + currentClassification+ 'andVirus_'+ filteredBy_EVEs + '_FilteredBy_'+ filteredBy_TEs + '_' + currentOrientation + 'Strand_heatmap_' + analysisType + '.pdf', dpi=600)

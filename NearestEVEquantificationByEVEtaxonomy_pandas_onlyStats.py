# -*- coding: utf-8 -*-
"""
@author: zwhitfield

Made to do the stats portion only for TEs near and/or overlapping EVEs
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:04:24 2016

@author: zwhitfield
This script takes output of NearestEVEquantification_pandas_overlapOrNearest_createFiles.py (need to use the Taxonomy
 files as input )and NearestEVEquantification_GenomeWide_pandasBash_NowWithStats_FrozenDataNoTEfam.py

Tests significance of enrichment of certain TE types near EVEs (derived from a particular viral taxonomy) compared to
their prevelance in the entire genome

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
plt.ion()
#matplotlib.__version__
plt.style.use('ggplot')

inputdir = str(sys.argv[1])
outputdir = str(sys.argv[2])

filteredBy = str(sys.argv[3])
filteredByCategory = str(sys.argv[4])
groupedCategory = str(sys.argv[5])
filteredBy_EVEs = str(sys.argv[6])

analysisType = str(sys.argv[7])

#inputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/"
#outputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/publicReleaseOutputCheck/forStats/"
#filteredBy = "NONE" #Enter 'NONE' if don't want any filtering.
#filteredByCategory = "TEclass"
#groupedCategory = "TEclass"
#filteredBy_EVEs = "family"

possibleOrientations = ['same','opposite']
currentClassification = groupedCategory

#will not be graphing, but want to be consistent with other scripts
distanceCutoffs = [-20000,20000]
print "Quantifying upstream and downstream TEs closest to EVEs"

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------FUNCTIONS------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def getHistData (NearestTEdata, orientation, EVEtax, typeOfAnalysis):
    if orientation == 'same':
        upANDdownstreamMERGED = NearestTEdata[NearestTEdata['EVEstrand'] == NearestTEdata['TEstrand']]
    if orientation == 'opposite':
        upANDdownstreamMERGED = NearestTEdata[NearestTEdata['EVEstrand'] != NearestTEdata['TEstrand']]

    #Get all taxonomy assignments under viral taxonomy level specified by EVEtax (ie filteredBy_EVEs).
    #so if EVEtax = 'family', will get list of EVEs all viral families in dataset.
    EVEclassifications = upANDdownstreamMERGED[EVEtax].unique()

    #loop through each instance of viral taxonomy for the given category. Stats will be determined for each type.
    #ie LTRs nearest EVEs derived from Flaviviridae. LTRs nearest EVEs derived from Rhabdoviridae etc...
    for currentEVEcategory in EVEclassifications:
        if (str(currentEVEcategory) == 'nan'):
            currentEVEcategory = 'Unclassified'
        upANDdownstreamMERGED_currentEVEcategory = upANDdownstreamMERGED[upANDdownstreamMERGED[EVEtax] == str(currentEVEcategory)]
        upANDdownstreamMERGED_currentEVEcategory = upANDdownstreamMERGED_currentEVEcategory.reset_index()

        print ("Filtering TEs by " + filteredBy + " in the " + orientation + " orientation as EVEs" + " and grouping by " + currentClassification)

        #will become list of arrays of distances (one per category of currentClassification)
        #distances will only be -1 0r 1 in this case (of only overlaps)
        distances=list()

        #will become list of names of each category of currentClassification
        names = list()

        #Use pandas groupby command to group data frame by shared currentClassification (ie group by entries with same 'TEdescription')
        groupedDF = upANDdownstreamMERGED_currentEVEcategory.groupby(currentClassification)

        #Obtain counts of each transposon in the entire genome, as produced by NearestEVEquantification_GenomeWide_pandas_NowWithStats.py
        wholeGenomeCounts = pd.read_csv(inputdir + 'ClassifiedBy_' + currentClassification+ '_FilteredBy_'+ filteredBy + '_WholeGenome.txt', sep='\t')

        #Get TOTAL number of TEs in genome
        wholeGenomeTotalCounts = sum(wholeGenomeCounts['Counts'])

        outfile = open(outputdir + 'ClassifiedBy_' + currentClassification + '_FilteredBy_' + filteredBy +
                       'AND' + currentEVEcategory + '_StatsComparedToWholeGenome_' + orientation +
                       'Strand_' + typeOfAnalysis + '.txt', 'w')

        outfile.write(currentClassification + '\t' +
                      "specificTEcounts" + '\t' +
                      "totalTEcountsOfType_" + filteredBy + '\t' +
                      "pValueOneSidedBinom" + '\t' +
                      "pValueOneSidedFEtest" + "\n")

        #Loop through groupedDF, getting the name (ie specific entry in currentClassification) and group (data in dataframe)
        #From there, extract distances from each group, append them to variable 'distance'. Do same for name
        #At the end, have 2 lists: each entry in 'distances' is all distances from EVEs of that currentCategory and each entry in 'names' is name of respective category in 'distances'
        #Example: If filteredBy was set to 'LTR', and classifications was set to 'TEfamily'
        #names=['Ty3_gypsy', 'Ty1_copia','Pao_Bel']
        #distances = [[all distances from EVEs of Ty3_gypsy elements], [all distances from EVEs of Ty1_copia elements], [all distances from EVEs of Pao Bel elements]]
        for name,group in groupedDF:
            distances.append(group['Distance'])
            names.append(name)
            #In order to perform enrichment tests, get number of TEs of given category ('names') in whole genome.
            specificTEcounts = wholeGenomeCounts[wholeGenomeCounts[currentClassification]==name].Counts.item()
            bTest = sp.binom_test(x=len(group['Distance']),n=len(upANDdownstreamMERGED_currentEVEcategory['Distance']),p=float(specificTEcounts)/float(wholeGenomeTotalCounts),alternative='greater')
            oddsratio, FEpValue = sp.fisher_exact([[len(group['Distance']), specificTEcounts],
                                       [len(upANDdownstreamMERGED_currentEVEcategory['Distance'])-len(group['Distance']), wholeGenomeTotalCounts-specificTEcounts]],
                                        alternative = 'greater')
            outfile.write(name + '\t' + str(len(group['Distance'])) + '\t' + str(len(upANDdownstreamMERGED_currentEVEcategory['Distance'])) + '\t' + str(bTest) + '\t' + str(FEpValue) + "\n")
        outfile.close()

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#------------------------------END_OF_FUNCTIONS-----------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#------------------------------Read in EVE-TE pairs-----------------------------
#-------------------------------------------------------------------------------

filePath = inputdir + "TEsClosestToEVEs_" + analysisType + "_withEVEtaxonomy.txt"
TEsNearestEVEs = pd.read_csv(filePath,
                             sep="\t")

#-------------------------------------------------------------------------------
#------------------------------Filter and calculate stats-----------------------
#-------------------------------------------------------------------------------

#Filter based on distance cutoffs (set at beginning of script)
TEsNearestEVEs = TEsNearestEVEs[(TEsNearestEVEs['Distance']>=distanceCutoffs[0]) & (TEsNearestEVEs['Distance']<=distanceCutoffs[1])]

#--------------------------------------------------------------------------------
#If desired, filter by desired class, establisehd by filteredBy variable at top of script.
if filteredBy!= 'NONE':
    TEsNearestEVEs = TEsNearestEVEs[TEsNearestEVEs[filteredByCategory]==filteredBy]

#Perform enrichment tests
getHistData(TEsNearestEVEs, 'same', filteredBy_EVEs, analysisType)
getHistData(TEsNearestEVEs, 'opposite', filteredBy_EVEs, analysisType)

#!/bin/bash

SCRIPTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/ScriptsFinalVersions/WithArguments/publicRelease/forHeatmapQuantificationAndGraphing"

INPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/" #REMEMBER TO CHANGE THE SECOND INSTANCE OF THIS (BELOW) IF NECESSARY

FILTEREDBY="Ty3_gypsy" #Enter 'NONE' if don't want any filtering. Specifies category to filter TEs by. Will only keep TEs of this type for analysis.
FILTERDBYCATEGORY="TEfamily" #So script knows what 'column' of data frame to filter by
GROUPBY="TEdescription" # Specifies what level/category to quantify TEs by (ie what goes into groupBy function of pandas) for statistics and plotting on resulting histogram. Use 'CombinedGroup' if grouping by combined categories. Need to specify various categories to be combined in actual python script.

EVETAXCATEGORY="family" #'family', 'species', 'genus', etc...Need to make sure the 'family' variable in NearestEVEbyTaxonomyBash_pandasHeatmap_FrozenData.py is properly set to reflect this. These are how virus groups will be organized on heatmap.

#Which type of EVE-TE pairings should be analyzed?
# "overlapOrNearest": If an EVE has an overlapping TE(s), then use those pairs. If no overlapping TE, then look at 'nearest neighbor' TE (both upstream AND downstream).
# "nearestOnly": Ignore all TEs which overlap EVEs, and find 'nearest neighbor' TEs for all EVEs, both upstream and downstream.
# "overlapOnly": Only look at EVEs with a TE whose coordinates overlap. If no overlapping TE, that EVE is not used in the analysis.
ANALYSISTYPE="overlapOrNearest" # "overlapOrNearest" OR "nearestOnly" OR "overlapOnly".

#Generate stats of TEs near EVEs (as compared to genome-wide prevelance of given TE category)
OUTPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/forStats/"
python ${SCRIPTDIRECTORY}/"NearestEVEquantificationByEVEtaxonomy_pandas_onlyStats.py" ${INPUTDIRECTORY} ${OUTPUTDIRECTORY} ${FILTEREDBY} ${FILTERDBYCATEGORY} ${GROUPBY} ${EVETAXCATEGORY} ${ANALYSISTYPE}


INPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/forStats/"
OUTPUTDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Figures/Heatmaps/"
python ${SCRIPTDIRECTORY}/"NearestEVEbyTaxonomyBash_pandas_Heatmap_FrozenData.py" ${INPUTDIRECTORY} ${OUTPUTDIRECTORY} ${FILTEREDBY} ${GROUPBY} ${EVETAXCATEGORY} ${ANALYSISTYPE}


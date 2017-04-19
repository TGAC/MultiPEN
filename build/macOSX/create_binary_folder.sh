#!/bin/sh
# Script to create folder for current version of MultiPEN
# It copies the following files:
# - MATLAB files for current binary
# - R scripts for pathway analysis
# - Scripts to test MultiPEN with example files
# To a new folder for MultiPEN's current version (i.e., MultiPEN_current-version_OS)

# Original files are located in:
#   Binary: "MultiPEN/MultiPEN_executable_development/mcc_files/"
#   R scripts for pathway analysis: "MultiPEN/scripts/pathwayAnalysis/"
#   Bash scripts to run examples: "MultiPEN/scripts/tests/"


# Currrent Version: MultiPEN_v003_OS
# Content:
#   MultiPEN (the application)
#   compileNetworkStringDB.R
#   pathwayAnalysis/
#      enrichmentGO_sortedList.R
#      enrichmentGO.R
#   example_cross_validation.sh
#   example_feature_selection.sh
#   example_pca.sh
#   example_hierarchical_clustering.sh
#   example_enrichment_GO.sh
#   example_enrichment_KEGG.sh
#   example_STRINGdb.sh

# folderName is the target folder for current version of MultiPEN
# change the following accordingly
folderName="MultiPEN_v003_OS"

origin="mcc_files/"   # directory containing binary
rscripts="../../scripts/pathwayAnalysis/"  # directory containing all R scripts for pathway analysis
exampleScripts="../../scripts/tests/OS/"  #path to bash scripts to run examples
stringdbRscript="../../scritps/compileNetworkStringDB.R"  # R script to compile network with STRINGdb

if [ ! -d "$folderName/" ]
then
   mkdir $folderName
else  # delete content
   rm -rf ${folderName}/MultiPEN*
   rm -r ${folderName}/*
fi

echo "###  Create Binary Folder  ###"
echo "Copying files from  ${origin}   to    ${folderName}/"

# Copy MultiPEN application
cp -Rp ${origin}MultiPEN.app MultiPEN_v001_OS/

# Create pathwayAnalysis folder
mkdir ${folderName}/pathwayAnalysis
# Copy all R scripts for pathway analysis
cp -Rp ${rscripts}*.R ${folderName}/pathwayAnalysis/  # copy all R scripts
# Copy R script to compile network using STRINGdb
cp -Rp $stringdbRscript ${folderName}


# Copy all scripts to run examples
cp ${exampleScripts}example_*.sh ${folderName}/  #copies all bash scripts to run examples 

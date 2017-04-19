#!/bin/sh
# Script to create folder for current version of MultiPEN
# It copies the following files:
# - MATLAB files for current binary
# - R scripts for pathway analysis
# - Scripts to test MultiPEN with example files
# To a new folder for MultiPEN's current version (i.e., MultiPEN_current-version_Linux)

# Original files are located in:
#   Binary: "MultiPEN/MultiPEN_executable_development/mcc_files_Linux/"
#   R scripts for pathway analysis: "MultiPEN/scripts/pathwayAnalysis/"
#   Bash scripts to run examples: "MultiPEN/scripts/tests/"


# Currrent Version: MultiPEN_v003_Linux
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
folderName="MultiPEN_v003_Linux"

origin="mcc_files_Linux/"   # directory containing binary
rscripts="../../scripts/pathwayAnalysis/"  # directory containing all R scripts for pathway analysis
exampleScripts="../../scripts/tests/Linux/"  #path to bash scripts to run examples
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
cp ${origin}MultiPEN ${folderName}
cp ${origin}run_MultiPEN.sh ${folderName}
chmod 755 ${folderName}/MultiPEN
chmod 755 ${folderName}/run_MultiPEN.sh

# Create pathwayAnalysis folder 
mkdir ${folderName}/pathwayAnalysis
# Copy all R scripts for pathway analysis
cp -Rp ${rscripts}*.R ${folderName}/pathwayAnalysis/  # copy all R scripts
# Copy R script to compile network using STRINGdb
cp -Rp $stringdbRscript ${folderName}


# Copy all scripts to run examples
cp -p ${exampleScripts}example_*.sh ${folderName}/  #copies all bash scripts to run examples
chmod 755 ${folderName}/example_*.sh

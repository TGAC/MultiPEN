#!/bin/sh
# Script to create folder for current version of MultiPEN
# It copies the following files:
# - MATLAB files for current binary
# - R scripts for pathway analysis
# - Scripts to test MultiPEN with example files
# To a new folder for MultiPEN's current version (i.e., MultiPEN_version_OS)

# Original files are located in:
#   Binary: "MultiPEN/MultiPEN_executable_development/mcc_files_version_OS/"
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

# Receive the operating system and version
if [ $# -eq 0 ] || [ $# -eq 1 ]
then
echo "Usage: $0 operating_system version";
exit 1;
else
OS=$1
version=$2
echo "Copying binary for $OS, version $version"
fi


# check the OS is supported
if [ "$OS" != "macOSX" ] && [ "$OS" != "Linux" ]; then
  echo "Currently, only supporting: macOSX and Linux"
  exit 1;
fi


# folderName is the target folder for current version of MultiPEN
# change the following accordingly
folderName="${OS}/MultiPEN_v${version}_${OS}"

origin="${OS}/mcc_files_v${version}_${OS}/"   # directory containing binary
rscripts="../scripts/pathwayAnalysis/"  # directory containing all R scripts for pathway analysis
exampleScripts="../scripts/tests/"  #path to bash scripts to run examples
stringdbRscript="../scripts/compileNetworkStringDB.R"  # R script to compile network with STRINGdb


if [ ! -d "$folderName/" ]
then
mkdir $folderName
echo "Creating folder for standalone application"
else  # delete content
rm -rf ${folderName}/MultiPEN*
rm -r ${folderName}/*
echo "Deleting old content... "
fi

echo "Creating Binary Folder ... "
echo "Copying files from  ${origin}   to    ${folderName}/"

# Copy MultiPEN application
if [ "$OS" == "macOSX" ]; then
   cp -Rp ${origin}MultiPEN.app ${folderName}
elif [ "$OS" == "Linux" ]; then
   cp ${origin}MultiPEN ${folderName}
   cp ${origin}run_MultiPEN.sh ${folderName}
   chmod 755 ${folderName}/MultiPEN
   chmod 755 ${folderName}/run_MultiPEN.sh
fi

# Create pathwayAnalysis folder 
mkdir ${folderName}/pathwayAnalysis
# Copy all R scripts for pathway analysis
cp -Rp ${rscripts}*.R ${folderName}/pathwayAnalysis/  # copy all R scripts
# Copy R script to compile network using STRINGdb
cp -Rp ${stringdbRscript} ${folderName}


# Copy all scripts to run examples
cp -p ${exampleScripts}example_*.sh ${folderName}/  #copies all bash scripts to run examples
if [ "$OS" == "Linux" ]; then
  chmod 755 ${folderName}/example_*.sh
fi

echo "Done!"

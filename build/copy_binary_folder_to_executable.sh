#!/bin/sh
# Script to copy MultiPEN's current version to MultiPEN_executable
# It copies the following files:
#   - MATLAB binary
#   - R scripts for pathway analysis
#   - Example files (input and output files) and
#   - Scripts to test MultiPEN with example files


# Source directory:
#   "MultiPEN/build/operating_system/MultiPEN_current-version_operating-system/"

# Target directory:
#   "MultiPEN/MultiPEN_executable/MultiPEN_current-version_operating-system/"


# Currrent Version: MultiPEN_v003_operating-system
# Content:
#   MultiPEN (the application)
#   pathwayAnalysis/
#      enrichmentGO_sortedList.R
#      enrichmentGO.R
#   compileNetworkStringDB.R
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

# folderName is the source folder for current version of MultiPEN
# change the following accordingly
folderName="/MultiPEN_v${version}_${OS}"


echo "###  Copy Binary to MultiPEN_Executable  ###"

# Copy folder for current version to MultiPEN_executable
target="../MultiPEN_executable/"
target=$target$folderName

echo "Copying files from  ${folderName}/   to  $target/"

# If folder exists, delete content, then copy updated content
if [ ! -d "$target" ]
then
mkdir $target
echo "${target}/ did not exist"
else  # delete content
rm -rf ${target}/MultiPEN*
rm -r ${target}/
echo "${target}/ and all its content will be been replaced with the updated version"
fi

cp -Rp ${folderName}/* ${target}/
echo "... done!"

# Compress folder for current version
echo 'Compressing the executable...'
zipFile="../MultiPEN_executable/MultiPEN_v${version}_${OS}.zip"
echo " ... into file: $zipFile"
zip -r -X $zipFile ${target}/
echo "Compressed file in: $zipFile"

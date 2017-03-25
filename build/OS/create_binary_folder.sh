#!/bin/sh
# Script to create folder for current version of MultiPEN
# It copies the following files:
# - MATLAB files for current binary
# - R scripts for pathway analysis
# - Example files (inputs and outputs) and
# - Scripts to test MultiPEN with example files
# To a new folder for MultiPEN's current version (i.e., MultiPEN_current-version_OS)

# Original files are located in:
#   Binary: "MultiPEN/MultiPEN_executable_development/mcc_files/"
#   R scripts for pathway analysis: "MultiPEN/scripts/pathwayAnalysis/"
#   Example Input Files: "MultiPEN/MultiPEN_executable_development/ExampleInputs/"
#   Example Output Files: "MultiPEN/MultiPEN_executable_development/ExampleOutputs/"
#   Bash scripts to run examples: "MultiPEN/scripts/tests/"


# Currrent Version: MultiPEN_v001_OS
# Content:
#   MultiPEN (the application)
#   pathwayAnalysis/
#      enrichmentGO_sortedList.R
#      enrichmentGO.R
#   example_cross_validation.sh
#   example_feature_selection.sh   
#   ExampleInputs
#      expressionData.txt
#      interactionMatrix.txt
#      sampleClass.txt
#   ExampleOutputs/
#      CrossValidation/
#         cross-validation_stats.txt
#      FeatureSelection/
#         MultiPEN-feature-selection_config.txt
#         MultiPEN-performance_feature-selection_lambda0.0001.txt
#         MultiPEN-Rankings_lambda0.0001_higher-in-cases.txt
#         MultiPEN-Rankings_lambda0.0001_higher-in-control.txt
#         MultiPEN-Rankings_lambda0.0001.txt
#         MultiPEN-vts_lambda0.0001.txt 

origin="mcc_files/"   # directory containing binary
rscripts="../../scripts/pathwayAnalysis/"  # directory containing all R scripts for pathway analysis
exampleScripts="../../scripts/tests/OS/"  #path to bash scripts to run examples


folderName="MultiPEN_v001_OS"   # target folder for current version of MultiPEN

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
cp -Rp ${rscripts}*.R ${folderName}/pathwayAnalysis/  # copies all R scripts

# Copy all example files and scripts to run examples
cp -Rp ExampleInputs/ ${folderName}/ExampleInputs/
cp -Rp ExampleOutputs/ ${folderName}/ExampleOutputs/
cp ${exampleScripts}example_*.sh ${folderName}/  #copies all bash scripts to run examples 

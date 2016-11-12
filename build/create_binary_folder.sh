#!/bin/sh
# Script to create folder for current version of MultiPEN
# It copies files for current binary, examples (input and output files) and 
# scripts to test the examples (original files) into a new folder for current version of MultiPEN (MultiPEN_current-version_OS)

# Original files are located in:
#   Binary: "MultiPEN/MultiPEN_executable_development/mcc_files/"
#   Example Input Files: "MultiPEN/MultiPEN_executable_development/ExampleInputs/"
#   Bash scripts to run examples: "MultiPEN/scripts/tests/"

# Currrent Version: MultiPEN_v001_OS
# Content:
#   MultiPEN (the application)
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

exampleScripts="../scripts/tests/"  #path to bash scripts to run examples
origin="mcc_files/"   # directory containing binary

folderName="MultiPEN_v001_OS"

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

# Copy all example files and scripts to run examples
cp -Rp ExampleInputs/ ${folderName}/ExampleInputs/
cp -Rp ExampleOutputs/ ${folderName}/ExampleOutputs/
cp ${exampleScripts}example_*.sh ${folderName}/  #copies all bash scripts to run examples 

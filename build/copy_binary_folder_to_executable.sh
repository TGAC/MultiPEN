#!/bin/sh
# Script to copy folder for current version to MultiPEN_executable
# It copies files for current binary, examples (input and output files) and 
# scripts to test the examples (original files) into a new folder for current version of MultiPEN (MultiPEN_current-version_OS)

# Original files are located in:
#   "MultiPEN/build/MultiPEN_current-version_OS/"


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


folderName="MultiPEN_v001_OS"

echo "###  Copy Binary to MultiPEN_Executable  ###"

# Copy folder for current version to MultiPEN_executable
target="../MultiPEN_executable/"
target=$target$folderName

echo "Copying files from  ./${folderName}/   to  $target/"

# If folder exists, delete content, then copy updated content
if [ ! -d "$target" ]
then
#mkdir $target
echo "${target}/ did not exist"
else  # delete content
rm -rf ${target}/MultiPEN*
rm -r ${target}/
echo "${target}/ and all its content has been replaced with new version"
fi

cp -Rp ${folderName} ${target}/

# Compress folder for current version
zipFile="../MultiPEN_executable/MultiPEN_v001_OS.zip"
zip -r -X $zipFile $folderName
echo "Compressed file in: $zipFile"

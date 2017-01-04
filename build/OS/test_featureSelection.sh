#!/bin/sh
# script to test the deployed MultiPEN application
# Module: Feature Selection

if [ $# -eq 0 ]
   then
      lambda=0.0001
   else
      lambda=$1
fi

MultiPEN="mcc_files/MultiPEN.app/Contents/MacOS/applauncher"
OutputDirectory="ExampleOutputs/"
ExpressionData="ExampleInputs/expressionData.txt"
Interactions="ExampleInputs/interactionMatrix.txt"
SampleClass="ExampleInputs/sampleClass.txt"

echo $MultiPEN
echo $lambda


# Run MultiPEN: feature selection 
# with default D (decision threshold) and number of iterations (numIter)

$MultiPEN FeatureSelection $OutputDirectory $ExpressionData $Interactions $SampleClass $lambda






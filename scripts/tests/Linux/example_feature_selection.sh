#!/bin/sh
# script to test the deployed MultiPEN application
# Module: Feature Selection

if [ $# -eq 0 ]
   then
      lambda=0.0001
   else
      lambda=$1
fi

MultiPEN="./run_MultiPEN.sh"
#Change the path to the compiler accordingly
#mcrPath="/usr/local/MATLAB/R2015b"
mcrPath="/usr/local/MATLAB/MATLAB_Runtime/v90"
OutputDirectory="ExampleOutputs/tests/"
ExpressionData="ExampleInputs/expressionData.txt"
Interactions="ExampleInputs/interactionMatrix.txt"
SampleClass="ExampleInputs/sampleClass.txt"

# Run MultiPEN: feature selection 
# with default decision threshold (D) and number of interactions (numIter)

$MultiPEN $mcrPath FeatureSelection $OutputDirectory $ExpressionData $Interactions $SampleClass $lambda


#!/bin/sh
# script to test the deployed MultiPEN application
# Module: Feature Selection

if [ $# -eq 0 ]
   then
      lambda=0.0001
   else
      lambda=$1
fi

MultiPEN="MultiPEN.app/Contents/MacOS/applauncher"
OutputDirectory="ExampleOutputs/"
ExpressionData="ExampleInputs/expressionData.txt"
Interactions="ExampleInputs/interactionMatrix.txt"
SampleClass="ExampleInputs/sampleClass.txt"
Folds=3

# Run MultiPEN: feature selection 
# with default number of interactions (numIter)

$MultiPEN CrossValidation $OutputDirectory $ExpressionData $Interactions $SampleClass $lambda $Folds


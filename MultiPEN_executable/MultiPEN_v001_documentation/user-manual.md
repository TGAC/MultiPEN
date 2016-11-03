# Running MultiPEN

## Cross Validation

### Example - OS

In the command line, navigate to the folder where the binary for MultiPEN is located, i.e., binary-OS/MultiPEN_v001_OS/. Then create variables for the paths to stand-alone application, output directory and input files by typing:

```
MultiPEN="MultiPEN.app/Contents/MacOS/applauncher"
OutputDirectory="ExampleOutputs/"
ExpressionData="ExampleInputs/X.txt"
Interactions="ExampleInputs/E.txt"
SampleClass="ExampleInputs/Y.txt"
lambda=0.001
Folds=3
NumIter=100
```

Note that in this example we are using the example files provided with the application. All the files used as example as located in the folder: ExampleInputs/.

Next, run Cross Validation with the following command:

```
$MultiPEN crossValidation $OutputDirectory $ExpressionData $Interactions $SampleClass $lambda $Folds $NumIter
```

To test more than one lambda one can specify the lambdas by using the notation (include the quotation mark symbols): 
“[lambda1 lambda2 … lambda3]”

For example, if we want to try two lambdas, namely 0.02 and 0.2, we would use the following command:

```
$MultiPEN crossValidation $OutputDirectory $ExpressionData $Interactions $SampleClass "[0.02 0.2]" $Folds $NumIter
```

# MultiPEN

MultiPEN uses a penalised logistic regression approach to find a subset of features (genes and/or metabolites) that hold more discriminant power to separate two classes, namely control and cases [1]. Such approach uses a molecular interaction network (e.g., protein-protein interaction network or metabolic network) to find the largest connected component that best separates the two conditions (for details on the logistic regression program to be optimised refer to [1]).

References

[1] GenePEN: analysis of network activity alterations in complex diseases via the pairwise elastic net., Vlassis N, Glaab E., Stat Appl Genet Mol Biol. 2015 Apr;14(2):221-4. doi: 10.1515/sagmb-2014-0045.


## Getting Started

MultiPEN is shared as a MATLAB stand-alone application, which requires the installation of the MATLAB Runtime for R2015b in your system. 

1.	Download and save MATLAB Runtime for R2015b for your operative system which can be found from:
http://www.mathworks.com/products/compiler/mcr/index.html 

2.	Double click the installer and follow the instructions in the installation wizard.




## Cross Validation

A common practice in the machine learning community is to first solve for the  parameter that optimises the logistic regression problem in equation 1 for your specific data. In MultiPEN, the module to precisely do that is crossValidation. 

### Syntax

*MultiPEN*  **crossValidation** *OutputDirectory ExpressionData Interactions SampleClass lambdas Folds NumIterations*


### Description


Parameter | Description
----------|-------------
**MultiPEN** | This is the path to the binary executable of MultiPEN, i.e., binary-OS/MultiPEN_v001_OS/.
**OutputDirectory** | Specify directory for output files.
**ExpressionData** |  The expression data is in tabular format where the rows are the features (genes and/or metabolites) and the columns are the samples. An example of a file containing expression data is shown in Figure 1 c).
**Interactions** |  The interaction matrix where the ith interaction (row) is represented as: [source target score] where *source* and *target* are names (symbolID for genes and CHEBI IDs for metabolites) of the connected nodes and *score* is a number in the range [0,1] representing the interaction confidence (where 1 corresponds to the maximum level of confidence). An example is shown in Figure 1 b).
**SampleClass** | For each sample specify if control (0) or case (1). An example of this file is shown in Figure 1 a) where each row contains the class for one sample. 
**lambdas** | Set of lambdas to test for cross validation. If wanting to test more than one lambda, specify the lambdas by using the notation (include the quotation mark symbols): “[lambda1 lambda2 … lambda3]”. For example, if we want to try two lambdas, namely 0.02 and 0.2, we would specify it with: “[0.02 0.2]”.
**Folds** | For cross validation.
**NumIterations** | Maximum number of iterations for the optimisation solver. Default value is 100.



![example inputs](images/figure-example-input-files.png)

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

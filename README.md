# MultiPEN

MultiPEN is a tool for the analysis of multi-omics data. We are developing our first version which will include a tool to analyse transcriptomics and metabolomics data. We are extending an approach based on a regularised logistic regression model [1] that uses a biological network as a priori knowledge.


### Using the Application

MultiPEN is shared as a MATLAB stand-alone application, which requires the installation of the MATLAB Runtime:

1. Download and save MATLAB Runtime for R2015b for your operating system [here](http://www.mathworks.com/products/compiler/mcr/index.html). 

2. Double click the installer and follow the instructions in the installation wizard.

3. Download the application (v.0.0.1 to be released soon! - OS and Linux).

#### Running MultiPEN

Details on how to use MultiPEN can be found in the [documentation](/MultiPEN_executable/MultiPEN_v003_documentation/user-manual.md).


### Using the MATLAB code

To use the MATLAB code, install the required libraries:

1. Install [TFCOS](http://cvxr.com/tfocs/) 
2. Install the MATLAB [gaimc package](http://www.mathworks.com/matlabcentral/fileexchange/24134-gaimc) (to determine connected network components).
3. Download [GenePEN](http://lcsb-portal.uni.lu/software/index.html) [1]. (Note that in this work we use a slightly modified version of GenePEN. More details to be provided soon). 


## About Us
This repository is developed by the [Jurkowski Group](http://www.earlham.ac.uk/jurkowski-group) in the Earlham Institute. For any questions or comments please contact Wiktor at wiktor.jurkowski@earlham.ac.uk. 



## References
1. GenePEN: analysis of network activity alterations in complex diseases via the pairwise elastic net., Vlassis N, Glaab E., Stat Appl Genet Mol Biol. 2015 Apr;14(2):221-4. doi: 10.1515/sagmb-2014-0045.


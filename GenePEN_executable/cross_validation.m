%function cross_validation(XDataFile, EDataFile, YDataFile, all_lambdas,
%cvfold, numIter, dirResults) deployed function
function [weights, vts, Fcross, outcome_stats] = ...
    cross_validation(XDataFile, EDataFile, YDataFile, all_lambdas, cvfold, numIter, dirResults)

% Example application of the GenePEN method for
% supervised sample classification and identification of
% network-grouped predictive features in gene/protein
% expression data
% 
% This script runs a cross-validation on example
% gene expression and molecular interaction data for
% GenePEN and computes a matrix of averaged outcome
% statistics, including the avg. size of the largest connected
% component and the avg. areq under the receiver operating
% characterisitc cure (AUC).

% Inputs:
%       XDataFile
%       EDataFile
%       YDataFile
%       all_lambdas     select lambda parameters to be used
%                       (here either a single lambda parameter can be
%                       specified or multiple lambda parameters can be
%                       tested) e.g. 
%                       all_lambdas = logspace(-12,2,num_of_lambdas);  
%                       all_lambdas = logspace(-0.3,-0.01,num_of_lambdas);
%       cvfold          the number of cross-validation cycles
%       numIter         maximum number of itrations for optimisation
%       dirResults      output directory
%       
% Outputs:

% Installation instructions: 
% Before running the code below, please add the path to
% your local TFOCS installation to Matlab and set the current
% directory to the one containing the gene/protein
% expression and interaction network data and file GenePEN.m.
% In order to compute the largest connected component
% among the selected features in the network, the Matlab
% package 'gaimc' has to be installed and the path to the
% must be added (gaimc is freely available online at
% www.mathworks.com/matlabcentral/fileexchange/24134-gaimc).
%
% Minimum requirements:
% The code has been tested on a standard laptop with 4 GB
% RAM and 2 GHz CPU. We recommend to use a machine with
% similar or superior memory and processor configuration.


% Please uncomment and edit the two lines below so that they
% point to your TFOCS and gaimc installation directories
%addpath('Libraries/GenePEN/'])
%addpath('Libraries/TFOCS-1.3.1/')
%addpath('Libraries/gaimc/')


% Read the example data
% X: data matrix with rows = samples, columns = genes/proteins
% E: interaction matrix with rows and column corresponding to genes/proteins
% Y: outcome label vector (1 for patient samples, 0 for controls)

%Input arguments passed from the system prompt will be received as strings
%  Thus, converting strings to double if required 
%X
if isa(XDataFile,'double')
    X = XDataFile;     
else
    X = load(XDataFile, '-ascii');
end
%E - network (edges)
if isa(EDataFile,'double')
    E = EDataFile;    
else
    E = load(EDataFile, '-ascii');
end
%class
if isa(YDataFile,'double')
    Y = YDataFile;
else
    Y = load(YDataFile, '-ascii');
end
%lambdas for cross validation
if ~isa(all_lambdas, 'double')
    all_lambdas = str2num(all_lambdas);   
end
%folds for cross validation
if ~isa(cvfold, 'double')
    cvfold = str2num(cvfold);
end
%maximum number of iterations for optimisation
if ~isa(numIter, 'double')
    numIter = str2num(numIter);
end


disp('lambdas for cross validation')
disp(all_lambdas)

disp('number of folds for cross validation: ')
disp(cvfold)

fprintf('Number of max iterations: \n\t%i\n\n', numIter)

disp('Output Directory: ')
disp(dirResults)


% store interaction matrix in sparse format
P = sparse(E(:,1),E(:,2),E(:,3),size(X,2),size(X,2),size(E,1));

% convert triangle matrix to symmetric matrix
A = P + P' - diag(diag(P));

% pre-process input matrix (For example Data)
%X = X - repmat(mean(X),size(X,1),1); 
%X = X / max(vec(X));

%
% Compute a cross-validation for GenePEN
F = zeros(length(all_lambdas),6);
for i = 1:length(all_lambdas)

    lambda = all_lambdas(i);
    fprintf('\n ######Lambda: %d\n', lambda)

    % use this line to ensure reproducibility
    % rand('seed',12345);  %PTR: it is an obsolote use as for MALTAB20015b
    rng('default')  %added by PTR
    
    
    if cvfold>1         
        c = cvpartition(Y,'kfold',cvfold);
    end   

    Fcross = zeros(cvfold,3);

    weights = zeros(numel(X(1,:)),cvfold);
    vts = zeros(1,cvfold);
    for j=1:cvfold,

            if cvfold==1        
                Xtrain = X;
                Ytrain = Y;
            else
                Xtrain = X(c.training(j),:);
                Ytrain = Y(c.training(j),:);
            end


            [wt, vt] = GenePEN( Xtrain, Ytrain, A, lambda, numIter );
            weights(:,j) = wt;   %added by PTR
            vts(1,j) = vt;       %added by PTR

            S = find(abs(wt)>1e-8);


            %% Feature Selection
            if (size(S, 1) > 0)
                [~,p] = largest_component(A(S,S));
            else
                p = 0;
            end

            if cvfold==1
                Xtest = X;
                y = Y;
            else
                Xtest = X(c.test(j),:);
                y = Y(c.test(j),:);
            end

            y_pred = Xtest *  wt;

            if (max(y_pred) - min(y_pred) == 0)
              if(max(y_pred) < 0)
                ypred = zeros(1, size(y_pred,1));
              else
                ypred = ones(1, size(y_pred,1));
              end
           else
              ypred = (y_pred - min(y_pred))/(max(y_pred) - min(y_pred));
           end

           r = tiedrank(ypred);
           auc = (sum(r(y==1)) - sum(y==1) * (sum(y==1)+1)/2) /( sum(y<1) * sum(y==1));

           Fcross(j,:) = [sum(p),numel(S),auc];

           disp(' ');
           disp('Current cycle:');
           disp(horzcat('Area under the curve (AUC): ',num2str(auc)));
           disp(horzcat('Size of largest connected component (LCC): ',num2str(sum(p))));
           disp(' ');       
           
           
    end

    % averaged performance statistics
    avgperf = mean(Fcross(:,3));
    sdperf = std(Fcross(:,3));
    avgclus = mean(Fcross(:,1));
    avgsel = mean(Fcross(:,2));
    sdclus = std(Fcross(:,1));

    disp(' ');
    disp('-------------');
    disp(horzcat('Average AUC for selected lambda: ',num2str(avgperf)));
    disp(' ');

    F(i,:) = [lambda,avgclus,sdclus,avgsel,avgperf,sdperf];
    
    %% Save feature feature selection results for current lambda
    %check if drectory exists
    if exist(dirResults, 'dir') ~= 7
        mkdir(dirResults)
    end    
    %save results
    fileName = [dirResults 'weights_lambda' num2str(lambda) '.mat'];
    fprintf('Saving results in file: \n    %s\n', fileName)    
    save(fileName, 'weights')
    dlmwrite([dirResults 'weights_lambda' num2str(lambda) '.txt'],weights,'delimiter','\t');
    %vts parameter
    fileName = [dirResults 'vts_lambda' num2str(lambda) '.mat'];
    fprintf('Saving results in file: \n    %s\n', fileName)    
    save(fileName, 'vts')
    dlmwrite([dirResults 'vts_lambda' num2str(lambda) '.txt'],vts,'delimiter','\t');

end

% the resulting matrix "outcome_stats"
% contains the following performance statistics:
outcome_stats = F;
disp(' ');
disp('   lambda     LCC     std-LCC-size	features   AUC     std-AUC')    
disp(outcome_stats)
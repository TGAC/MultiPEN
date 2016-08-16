function [weights, vts, Fcross, outcome_stats] = ...
    cross_validation(X, E, Y, all_lambdas, cvfold, numIter)

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
%   XDataFile       data matrix with rows = samples, columns = genes/proteins
%                   assumes that data has been normalised/standardised
%   YDataFile       class for each sample, i.e., 0 for control, 1 for cases
%   EDataFile       Interaction matrix, i.e., the network's edge list. 
%                   Each row specifies an interaction in the format:
%                   [node node interaction-weight]
%   all_lambdas     select lambda parameters to be used
%                   (here either a single lambda parameter can be
%                   specified or multiple lambda parameters can be
%                   tested) e.g.:
%                   all_lambdas = logspace(-12,2,num_of_lambdas);  
%                   all_lambdas = logspace(-0.3,-0.01,num_of_lambdas);
%       cvfold      the number of cross-validation cycles
%       numIter     maximum number of itrations for optimisation
%       
% Outputs:
%   weights
%   vts  
%   Fcross
%   outcome_stats
%

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

%% Show input parameters in screen
disp('lambdas for cross validation')
disp(all_lambdas)

disp('number of folds for cross validation: ')
disp(cvfold)

fprintf('Number of max iterations: \n\t%i\n\n', numIter)


%% Preparing Data
% store interaction matrix in sparse format
P = sparse(E(:,1),E(:,2),E(:,3),size(X,2),size(X,2),size(E,1));

% convert triangle matrix to symmetric matrix
A = P + P' - diag(diag(P));


%% Cross Validation
% Compute a cross-validation for GenePEN
F = zeros(length(all_lambdas),6);
% weights has m columns, where m = length(all_lambdas) * cvfold;
weights = zeros(numel(X(1,:)), length(all_lambdas)*cvfold );
vts = zeros(1, length(all_lambdas)*cvfold );

for i = 1:length(all_lambdas)
    lambda = all_lambdas(i);
    fprintf('\n ######Lambda: %d\n', lambda)

    % use this line to ensure reproducibility
    % rand('seed',12345);  %PTR: it is an obsolote use as for MALTAB20015b
    rng('default')  %added by PTR
    
    % If feature selection, i.e., more than one fold
    if cvfold>1         
        c = cvpartition(Y,'kfold',cvfold);
    end   

    Fcross = zeros(cvfold,3);
    
    for j=1:cvfold,

            if cvfold==1        
                Xtrain = X;
                Ytrain = Y;
            else
                Xtrain = X(c.training(j),:);
                Ytrain = Y(c.training(j),:);
            end


            [wt, vt] = GenePEN( Xtrain, Ytrain, A, lambda, numIter ); 
            
            % indexes of features that would be selected
            S = find(abs(wt)>1e-8);
            
            % Results for GenePEN
            %index for the current column considering folds and lambdas
            indx = j + ((i-1)*cvfold) ;  %added by PTR
            weights(:,indx) = wt;   %added by PTR
            vts(1,j) = vt;       %added by PTR            
            weights(abs(wt)<1e-8,indx) = 0;

            
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

end

% the resulting matrix "outcome_stats"
% contains the following performance statistics:
outcome_stats = F;
disp(' ');
disp('   lambda     LCC     std-LCC-size	features   AUC     std-AUC')    
disp(outcome_stats)
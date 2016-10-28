function [weights, vt] = runGenePEN(X, E, Y, lambda, numIter)
% This function runs GenePEN [1] on expression data (gene expression 
% and/or metabolite level data) and a molecular interaction matrix.
% 
% [1] Vlassis N and Glaab E., Stat. Appl. Genet Mol Biol. 2015, 
%     doi: 10.1515/sagmb-2014-0045)

% Inputs:
%   X           expression matrix with rows = samples, columns = features
%               features are genes and/or metabolites
%               We assume that data has been normalised/standardised
%   E           Interaction matrix, i.e., the network's edge list. 
%               Each row specifies an interaction in the format:
%               [nodeA nodeB interaction-weight]
%   Y           sample class for each sample: 0 for control, 1 for cases
%   lambda      tradeoff parameter for the regularisation function
%   numIter     maximum number of iterations for optimisation
%       
% Outputs:
%   weights     weights calculated by GenePEN for each fearture
%   vts         the intercept
%   Ypred       predicted classes for samples in test set

%% Preparing Data
% store interaction matrix in sparse format
P = sparse(E(:,1),E(:,2),E(:,3),size(X,2),size(X,2),size(E,1));

% convert triangle matrix to symmetric matrix
A = P + P' - diag(diag(P));

[weights, vt] = GenePEN( X, Y, A, lambda, numIter ); 
weights(abs(weights)<1e-8) = 0;

function [FS, vts] = featureSelection(X, E, Y, featureAnnot, lambda, numIter)
% Feature selection using MultiPEN


% This function is called from: MultiPEN.m
% Outputs:
%           featureSelection.csv
%           which is a table with columns:
%           [name, weight, ranking]
% Inputs:
% X
% E
% Y
% lambda 
% outputDir
% numIter 

%% Feature Selection with MultiPEN 
% (number of folds is set to 1 to use all samples)
%[weights, ~, ~, ~] = cross_validation(X, E, Y, lambda, folds, numIter, outputDir);
[weights, vts, ~, ~] = cross_validation(X, E, Y, lambda, 1, numIter);
weights(weights < 1e-8) = 0;

%% Rank the feature selection
FS = table();
FS.name = featureAnnot;
FS.weight = weights;
R = tiedrank(abs(weights));
% Reverse the ranking order used in tiedrank
% so that the maximum absolute weight has ranking 1
% and the minimum absoulute weight has the larges ranking
n = numel(R);
FS.ranking = n - (R - 1);

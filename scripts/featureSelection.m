function FS = featureSelection(X, E, Y, Annot, lambda, outputDir, numIter)
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
[weights, ~, ~, ~] = cross_validation(X, E, Y, lambda, 1, numIter, outputDir);
weights(weights < 1e-8) = 0;

%% Rank the feature selection
FS = table();
FS.name = Annot;
FS.weight = weights;
R = tiedrank(abs(weights));
% Reverse the ranking order used in tiedrank
% so that the maximum absolute weight has ranking 1
% and the minimum absoulute weight has the larges ranking
n = numel(R);
FS.ranking = n - (R - 1);



%% Write feature selection (FS table) to file
%fileName = [outputDir 'Rankings_lambda' num2str(lambda)];
%fprintf('Writing feature selection to file: \n\t%s\n',fileName)
%writetable(FS, [fileName '.txt'], 'delimiter', '\t');
%save([fileName '.mat'], 'FS')
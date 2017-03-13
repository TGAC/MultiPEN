function saveFeatureSelection(outputDir, FS, vt, higherControl, higherCases, ...
    stats, edges, lambda, decisionThr, numIter, otherPar)
% WRITES feature selection (FS table) to file
% otherPar  is a string that contains information on selected preprocessing 
%           example: 010.500100
                    % Position      Parameter       Values
                    %   1         logTransform      0 or 1
                    %   2           normalise       0 or 1
                    %  3:6         decisionThr     1.00, 0.60, etc
                    %  7:10          numIter       0300, 1000, 2000, etc

%check if output directory exists
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir)
end


%feature's ranking
%[name weight ranking]
FileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) ...
    '_' otherPar '.txt'];
fprintf('Writing feature selection to file: \n\t%s\n',FileName)
writetable(FS, FileName, 'delimiter', '\t');

%intercept learnt from feature selection
FileName = [outputDir 'MultiPEN-vts_lambda' num2str(lambda) ...
    '_' otherPar '.txt'];
dlmwrite(FileName, vt, 'delimiter', '\t');    

%% write separate tables for higherControl an higherCases
FileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) ...
    '_higher-in-control_' otherPar '.txt'];
fprintf('Writing feature selection to file: \n\t%s\n',FileName)
writetable(higherControl, FileName, 'delimiter', '\t');

FileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) ...
    '_higher-in-cases_' otherPar '.txt'];
fprintf('Writing feature selection to file: \n\t%s\n',FileName)
writetable(higherCases, FileName, 'delimiter', '\t');

%% write table with statistics for feature selection accuracy
FileName = [outputDir 'MultiPEN-performance_feature-selection_lambda' ...
    num2str(lambda) '_' otherPar '.txt'];
fprintf('Writing performance for feature selection to file: \n\t%s\n',FileName)
writetable(stats, FileName, 'delimiter', '\t');

%% write file with the parameters used to run feature selection
config = table();
config.lambda = lambda;
config.numIter = numIter;
config.decisionThreshold = decisionThr;
FileName = [outputDir 'MultiPEN-feature-selection_config_lambda' ...
    num2str(lambda) '_' otherPar '.txt'];
writetable(config, FileName, 'delimiter', '\t');

%% Write network used for Feature Selection
FileName = [outputDir 'edges-for-feature-selection_lambda' ...
    num2str(lambda) '_' otherPar '.txt'];
writetable(edges, FileName, 'delimiter', '\t')
fprintf('\tEdges for feature selection are saved to file: \n\t%s\n', ...
    FileName)

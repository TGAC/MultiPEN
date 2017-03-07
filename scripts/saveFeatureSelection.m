function saveFeatureSelection(outputDir, FS, vt, higherControl, higherCases, ...
    stats, edges, lambda, decisionThr, numIter, logTransform, normalise, otherPar)

% WRITE feature selection (FS table) to file
% if ~strcmp(saveResults,'false')
%     if strcmp(saveResults, 'true')
%         outputDir = 'output_MultiPEN/feature_selection/';
%     else
%         outputDir = saveResults;
%     end

    %check if output directory exists
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir)
    end
    
    % include details on optional parameters at the end of file name
    optional = '';
    if logTransform && normalised
        optional = '_log2_normalised';
    elseif logTransform && ~normalise
        optional = '_log2';
    elseif ~logTransform && normalise
        optional = '_normalised';
    end
    
    %feature's ranking
    %[name weight ranking]
    FileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) ...
        optional '_' otherPar '.txt'];
    fprintf('Writing feature selection to file: \n\t%s\n',FileName)
    writetable(FS, FileName, 'delimiter', '\t');

    %intercept learnt from feature selection
    FileName = [outputDir 'MultiPEN-vts_lambda' num2str(lambda) ...
        optional '_' otherPar '.txt'];
    dlmwrite(FileName, vt, 'delimiter', '\t');    

    %% write separate tables for higherControl an higherCases
    FileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) ...
        '_higher-in-control' optional '_' otherPar '.txt'];
    fprintf('Writing feature selection to file: \n\t%s\n',FileName)
    writetable(higherControl, FileName, 'delimiter', '\t');

    FileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) ...
        '_higher-in-cases' optional '_' otherPar '.txt'];
    fprintf('Writing feature selection to file: \n\t%s\n',FileName)
    writetable(higherCases, FileName, 'delimiter', '\t');

    %% write table with statistics for feature selection accuracy
    FileName = [outputDir 'MultiPEN-performance_feature-selection_lambda' ...
        num2str(lambda) optional '_' otherPar '.txt'];
    fprintf('Writing performance for feature selection to file: \n\t%s\n',FileName)
    writetable(stats, FileName, 'delimiter', '\t');

    %% write file with the parameters used to run feature selection
    config = table();
    config.lambda = lambda;
    config.numIter = numIter;
    config.decisionThreshold = decisionThr;
    FileName = [outputDir 'MultiPEN-feature-selection_config_lambda' ...
        num2str(lambda) optional '_' otherPar '.txt'];
    writetable(config, FileName, 'delimiter', '\t');

    %% Write network used for Feature Selection
    FileName = [outputDir 'edges-for-feature-selection_lambda' ...
        num2str(lambda) optional '_' otherPar '.txt'];
    writetable(edges, FileName, 'delimiter', '\t')
    fprintf('\tEdges for feature selection are saved to file: \n\t%s\n', ...
        FileName)
% end

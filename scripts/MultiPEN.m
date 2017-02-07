function MP = MultiPEN(analysisType, saveResults, varargin)

% MultiPEN performs analysis of transcritpomics and metabolomics for feature selection
% Coypyright (C) 2016 {Perla Rey}
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% MultiPEN v0.0.1 uses:
% - GenePEN (Vlassis N and Glaab E., Stat. Appl. Genet Mol Biol. 2015, doi: 10.1515/sagmb-2014-0045)   
% - TFOCS-1.3.1
% - gaimc
% 


% INPUTS
%   analysisType    Options are: 
%                   'hierarchicalClustering'
%                   'CrossValidation', 'FeatureSelection'
%                   coming soon: 'enrichmentGO', 'RandomiseNetwork', 'ErdosRenyi'
%   ED              expression data with feature names and expression per
%                   sample, it's a p-by-(n+1) table 
%                   columns are: name, sample1, sample2, ..., sampleN
%   saveResults     Possible values are:
%                   'true' 'false' the_output_directory


%% VERIFY INPUT ARGUMENTS
switch analysisType
    case 'PCA'
        % PCA needs parameters:
        % expData (expression, samples, features), group, threshold (opt),
        % plotTitle (optional)  
        if (isempty(varargin)) || (length(varargin) < 2) || (length(varargin) > 4)
            error('The number of arguments is incorrect')
        else
            expData = varargin{1};            
            groups = varargin{2};
            switch length(varargin)
                case 3
                    threshold = str2num(varargin{3});
                    % plotTitle is optional and it won't be provided to HC function
                case 4
                    threshold = str2num(varargin{3});
                    plotTitle = varargin{4};
            end
        end
        
    case 'HierarchicalClustering'
        % Hierarchical clustering needs parameters:
        % D (expression, samples, features), saveFigure, threshold (opt),
        % plotTitle (optional)  
        if (isempty(varargin)) || (length(varargin) > 3)
            error('The number of arguments is incorrect')
        else
            expData = varargin{1};            
            switch length(varargin)
                case 2
                    threshold = str2num(varargin{2});
                    % plotTitle is optional and it won't be provided to HC function
                case 3
                    threshold = str2num(varargin{2});
                    plotTitle = varargin{3};               
            end
        end
        
    case 'CrossValidation'
        % cross validation needs parameters: 
        % D, E, Y, lambdas, folds, numIter (optional)
        if ~((length(varargin) == 5) || (length(varargin) == 6))
            error('The number of arguments is incorrect')
        else
            expData = varargin{1};   %expression data
            interactionMatrix = varargin{2};
            sampleClass = varargin{3};
            lambdas = str2num(varargin{4});
            folds = str2num(varargin{5});
            if length(varargin) == 5
                numIter = 100;
            else
                numIter = str2num(varargin{6});
            end
        end
        
    case 'FeatureSelection'
        % FeatureSelection needs parameters:
        % D, E, Y, lambda, decisionThr, numIter (optional)
        if ~((length(varargin) == 4) || (length(varargin) == 5) || ...
                (length(varargin) == 6))
            error('The number of arguments is incorrect')
        else
            expData = varargin{1};
            interactionMatrix = varargin{2};
            sampleClass = varargin{3};
            lambda = str2double(varargin{4});
            switch length(varargin)
                case 4  % values by default
                    decisionThr = 0.50;  %by default decision threshold
                    numIter = 100;    %by default number of iterations
                case 5
                    decisionThr = str2double(varargin{5});
                    numIter = 100;
                case 6
                    decisionThr = str2double(varargin{5});
                    numIter = str2num(varargin{6});
            end
            
        end
        
    case 'EnrichmentGO'
        % enrichmentGO needs parameters:
        % mpRankings (output from FeatureSelection: MultiPEN-Rankings_lambda{lambda}.txt)        
        if ~(length(varargin) == 1)
            error('The number of arguments is incorrect')
        else
            mpRankings = varargin{1};
        end
        
    otherwise
        error('Please specify a valid analysis')
end
    
%% READ INPUT DATA
%  Input arguments passed from the system prompt will be received as strings
%  Thus, converting strings to double if required 

% expression data is a table where:
%    the rows are the features (genes and/or metabolites)
%    the columns are the samples
% expression data can be provided as file or as a table
if exist('expData', 'var')
    if ~istable(expData)
        expData = readtable(expData, 'delimiter', '\t');
    end
    expression = expData;
    XAnnotation = expression.name;  % p-by-1  -  the features
    X = expression;
    X.name = [];   
    
    samples = X.Properties.VariableNames';  % n-by-1 
    X = table2array(X)';  % n-by-p
end

% E - interaction matrix (network's edges)
%     matrix can be provided as file or as a table
if exist('interactionMatrix', 'var')
    if ~istable(interactionMatrix)
        interactionMatrix = readtable(interactionMatrix, 'delimiter', '\t');
    end
end
    

% Sample Class
if exist('sampleClass', 'var')
    if ~istable(sampleClass)        
        sampleClass = readtable(sampleClass, 'delimiter', '\t');
    end
    Y = sampleClass.class;
end

% Groups for PCA
if exist('groups', 'var')
    if ~iscell(groups) 
        groups = readtable(groups, 'delimiter', '\t');
        groups = table2cell(groups)';
    end
end


%% ADD PATH TO LIBRARIES
if ~isdeployed
    addpath('Libraries/')
    addpath('Libraries/fastGapFill/')
    addpath('Libraries/gaimc/')
    addpath('Libraries/GenePEN/')
    addpath('Libraries/TFOCS-1.3.1/')
    addpath('Libraries/gscatter3/')
end


%% Analysis

switch analysisType
    case 'PCA'
        % X is the n-by-p data matrix for n samples, p features
        pcaExpressionData(X, groups, plotTitle, XAnnotation);
        
        MP = 1;  % exit code 1
    
    case 'HierarchicalClustering'
        %hierarchicalClustering(expression, samples, features, saveFigure, varargin)
        if ~ (exist('threshold','var') && exist('plotTitle', 'var') )
            hierarchicalClustering(X, samples, XAnnotation, saveResults)
        elseif (exist('threshold','var') && exist('plotTitle', 'var') )
            hierarchicalClustering(X, samples, XAnnotation, saveResults, threshold, plotTitle)
        elseif exist('threshold', 'var') 
            hierarchicalClustering(X, samples, XAnnotation, saveResults, threshold)
        elseif exist('plotTitle', 'var')
            hierarchicalClustering(X, samples, XAnnotation, saveResults, plotTitle)
        end
        
        MP = 1;   % exit code 1  
    
    case 'CrossValidation'                
        % Get subnetwork for  expressionData from interactionMatrix
        % i.e. use only interactions whose nodes correspond to features
        % in the expression data
        fprintf('##############\n')
        fprintf('Obtaining subnetwork for the expression data ... \n')        
        % edges has the edges for the subnetwork as table
        % with name of source, name of target, and weight
        [E, edges] = subnetwork4ExpressionData(interactionMatrix, XAnnotation);
        
        
        %CrossValidation for different lambdas
        fprintf('##############\n')
        fprintf('Performing cross validation... \n')        
        [~, ~,stats, yTest, yTestPred] = crossValidation(X, E, Y, lambdas, folds, numIter);
        
        if ~strcmp(saveResults,'false')
            if strcmp(saveResults, 'true')
                outputDir = 'output_MultiPEN/CrossValidation/';
            else
                outputDir = saveResults;
            end
            
            %check if output directory exists
            if exist(outputDir, 'dir') ~= 7
                mkdir(outputDir)
            end
                        
            % Statistics for cross validation
            fileName = [outputDir 'cross-validation_stats.txt'];
            writetable(stats, fileName, 'delimiter', '\t')
            fprintf('\tStatistics are saved to file: \n\t%s\n', ...
                fileName)
            fprintf('Following are the statistics for cross validation\n')
            display(stats)
                        
        end
        
        MP = stats;
        
          
    case 'FeatureSelection'
        % Get subnetwork for  expressionData from interactionMatrix
        % i.e. use only interactions whose nodes correspond to features
        % in the expression data
        fprintf('##############\n')
        fprintf('Obtaining subnetwork for the expression data ... \n')        
        E = subnetwork4ExpressionData(interactionMatrix, XAnnotation);
        
        %Feature selection for a specific lambda
        fprintf('Performing feature selection... \n')
        %FS is a table with columns: [name, weight, ranking]
        [FS, vt, stats] = featureSelection(X, E, Y, XAnnotation, lambda, numIter, decisionThr);
        
        %compute fold change
        [FS, higherControl, higherCases] = foldChange(FS, X, Y, samples);
        
        %sort results by ranking
        FS = sortrows(FS,{'ranking'},{'ascend'});
        
        % WRITE feature selection (FS table) to file
        if ~strcmp(saveResults,'false')
            if strcmp(saveResults, 'true')
                outputDir = 'output_MultiPEN/feature_selection/';
            else
                outputDir = saveResults;
            end
            
            %check if output directory exists
            if exist(outputDir, 'dir') ~= 7
                mkdir(outputDir)
            end
            
            %feature's ranking
            %[name weight ranking]
            fileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda)];
            fprintf('Writing feature selection to file: \n\t%s\n',fileName)
            writetable(FS, [fileName '.txt'], 'delimiter', '\t');
            
            %intercept learnt from feature selection
            fileName = [outputDir 'MultiPEN-vts_lambda' num2str(lambda)];            
            dlmwrite([fileName '.txt'], vt, 'delimiter', '\t');    
            
            %% write separate tables for higherControl an higherCases
            fileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) '_higher-in-control.txt'];
            fprintf('Writing feature selection to file: \n\t%s\n',fileName)
            writetable(higherControl, fileName, 'delimiter', '\t');
            
            fileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) '_higher-in-cases.txt'];
            fprintf('Writing feature selection to file: \n\t%s\n',fileName)
            writetable(higherCases, fileName, 'delimiter', '\t');
            
            %% write table with statistics for feature selection accuracy
            fileName = [outputDir 'MultiPEN-performance_feature-selection_lambda' num2str(lambda) '.txt'];
            fprintf('Writing performance for feature selection to file: \n\t%s\n',fileName)
            writetable(stats, fileName, 'delimiter', '\t');
            
            %% write file with the parameters used to run feature selection
            config = table();
            config.lambda = lambda;
            config.numIter = numIter;
            config.decisionThreshold = decisionThr;
            fileName = [outputDir 'MultiPEN-feature-selection_config.txt'];
            writetable(config, fileName, 'delimiter', '\t');
            
        end
                       
        MP = FS;
        
        
    case 'EnrichmentGO'
        % output directory 
        %if ~strcmp(saveResults, 'false')
            if strcmp(saveResults, 'true')
                outputDir = 'output_MultiPEN/enrichment-GO/';
            else
                outputDir = saveResults;
            end

            %check if output directory exists
            if exist(outputDir, 'dir') ~= 7 && ~strcmp(saveResults, 'false')
                mkdir(outputDir)
            end
        %end
        
        % Build string to call the R string enrichmentGO.R  (using Rscript)
        % syntaxis:
        % path_to_Rscript script_to_run file_name output_directory
        callToRscript = '/Library/Frameworks/R.framework/Resources/Rscript scripts/pathwayAnalysis/enrichmentGO.R';
        callToRscript = [callToRscript ' ' mpRankings ' ' outputDir];
        system(callToRscript)
        
        % Load the table with results
        MP = readtable([outputDir 'enrichment-GO.txt'], 'Delimiter', '\t');
        

end
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
        % D, E, Y, lambda, logTransform (optional), normalise (optional), 
        %      decisionThr (optional), numIter (optional)
        if ~((length(varargin) == 4) || (length(varargin) == 5))
            error('The number of arguments is incorrect')
        else
            expData = varargin{1};
            interactionMatrix = varargin{2};
            sampleClass = varargin{3};
            lambda = str2double(varargin{4});
            
            switch length(varargin)
                case 4  % values by default
                    logTransform = 0;
                    normalise = 0;
                    decisionThr = 0.50;  %by default decision threshold
                    numIter = 100;    %by default number of iterations
                case 5
                    % optional is a string, for example: 010.500100
                    % where digits of groups of digit represent the
                    % following parameter:
                    % (digit's) Position      Parameter       Values
                    %         1             logTransform      0 or 1
                    %         2               normalise       0 or 1
                    %        3:6             decisionThr     1.00, 0.60, etc
                    %        7:10              numIter       0300, 1000, 2000, etc
                    optional = varargin{5};
                    logTransform = str2num(optional(1));
                    normalise = str2num(optional(2));
                    decisionThr = str2double(optional(3:6));
                    numIter = str2double(optional(7:10));

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
        
        % Default set of 20 lambdas
        % If none were provided, i.e., lambdas=-1
        if lambdas == -1
            lambdas = logspace(-12,2,20);
        end
        
        % edges has the edges for the subnetwork as table
        % with name of source, name of target, and weight
        [E, edges] = subnetwork4ExpressionData(interactionMatrix, XAnnotation);
        [Xt, XAnnotationT] = genesInSubnetwork(X, XAnnotation, edges);
        [E, edges] = subnetwork4ExpressionData(edges, XAnnotationT);
        
        %CrossValidation for different lambdas
        fprintf('##############\n')
        fprintf('Performing cross validation... \n')        
        [~, ~,stats, yTest, yTestPred] = crossValidation(Xt, E, Y, lambdas, folds, numIter);
        
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
            FileName = [outputDir 'cross-validation_stats.txt'];
            writetable(stats, FileName, 'delimiter', '\t')
            fprintf('\tStatistics are saved to file: \n\t%s\n', ...
                FileName)
            fprintf('Following are the statistics for cross validation\n')
            display(stats)
            
            % Write network used for Cross Validation
            FileName = [outputDir 'edges-for-cross-validation.txt'];
            writetable(edges, FileName, 'delimiter', '\t')
            fprintf('\tEdges for cross validation are saved to file: \n\t%s\n', ...
                FileName)
                        
        end
        
        MP = stats;
        
          
    case 'FeatureSelection'
        % Get subnetwork for  expressionData from interactionMatrix
        % i.e. use only interactions whose nodes correspond to features
        % in the expression data
        fprintf('##############\n')
        fprintf('Obtaining subnetwork for the expression data ... \n')        
        [~, edges] = subnetwork4ExpressionData(interactionMatrix, XAnnotation);
        [Xt, XAnnotationT] = genesInSubnetwork(X, XAnnotation, edges);
        [E, edges] = subnetwork4ExpressionData(edges, XAnnotationT);
        
        % Remove duplicated edges
        E = removeDuplicatedEdges(E);
        
        %% Feature selection for a specific lambda
        fprintf('Performing feature selection... \n')
        % Data pre-processing
        if logTransform && normalise
            XtFS = zscore(log2(Xt));
        elseif logTransform && ~normalise
            XtFS = log2(Xt);
        elseif ~logTransform && normalise
            XtFS = zscore(Xt);
        elseif ~logTransform && ~normalise
            XtFS = Xt;
        end
        %FS is a table with columns: [name, weight, ranking]
        [FS, vt, stats] = featureSelection(XtFS, E, Y, XAnnotationT, lambda, numIter, decisionThr);
        
        %compute fold change
        [FS, higherControl, higherCases] = foldChange(FS, Xt, Y, samples);
        
        % add the expression data used for Feature Selection
        FS(:,end+1:end+(numel(XtFS(:,1)))) = array2table(XtFS');
        FS.Properties.VariableNames(6:end) = samples;
        
        %sort results by ranking
        FS = sortrows(FS,{'ranking'},{'ascend'});
        
        % how many features were selected
        fprintf('Number of features selected: %i\n', numel(FS.weight(FS.weight ~= 0)))
        
        %% WRITE feature selection (FS table) to file
        if ~strcmp(saveResults,'false')
            if strcmp(saveResults, 'true')
                outputDir = 'output_MultiPEN/feature_selection/';
            else
                outputDir = saveResults;
            end
            
            saveFeatureSelection(outputDir, FS, vt, higherControl, higherCases, ...
                stats, edges, lambda, decisionThr, numIter, optional);

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
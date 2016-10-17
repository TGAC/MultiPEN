function MP = MultiPEN(analysisType, saveResults, varargin)

% Function to perform analysis of omics data using MultiPEN and
% the different type of analysis (specified by parameter 'analysisType')
% MultiPEN 0.0.1 computes feature selection from transcritpomics and 
% metabolomics data. 
% Written by Perla Rey, Agust 2016

% MultiPEN extends the approach proposed by GenePEN (Vlassis N and Glaab E., 
% Stat. Appl. Genet Mol Biol. 2015, doi: 10.1515/sagmb-2014-0045)
%
% It uses following libraries:
% GenePEN
% TFOCS-1.3.1
% gaimc
% 

% INPUTS
%   analysisType    Options are: 
%                   'hierarchicalClustering'
%                   'crossValidation', 'featureSelection'
%                   coming soon: 'enrichmentGO', 'RandomiseNetwork', 'ErdosRenyi'
%   saveResults     Possible values are:
%                   'true' 'false' outputDirectory


%% VERIFY INPUT ARGUMENTS
switch analysisType
    case 'hierarchicalClustering'
        
        
    case 'crossValidation'
        % cross validation needs parameters: 
        % X, E, Y, lambdas, folds, numIter (optional)
        if ~((length(varargin) == 5) || (length(varargin) == 6))
            error('The number of arguments is incorrect')
        else
            X = varargin{1};
            E = varargin{2};
            Y = varargin{3};
            lambdas = str2num(varargin{4});
            folds = str2num(varargin{5});
            if length(varargin) == 5
                numIter = 100;
            else
                numIter = str2num(varargin{6});
            end
        end
        
    case 'featureSelection'
        % featureSelection needs parameters:
        %X, E, Y, lambda, features, sampleAnnot, numIter (optional)
        if ~((length(varargin) == 6) || (length(varargin) == 7))
            error('The number of arguments is incorrect')
        else
            X = varargin{1};
            E = varargin{2};
            Y = varargin{3};
            lambda = str2num(varargin{4});
            features = varargin{5};
            samples = varargin{6};
            if length(varargin) == 6
                numIter = 100;
            else
                numIter = str2num(varargin{7});
            end
        end
        
    case 'enrichmentGO'
        % enrichmentGO needs parameters:
        % fileName (output from featureSelection: MultiPEN-Rankings_lambda{lambda}.txt)        
        if ~(length(varargin) == 1)
            error('The number of arguments is incorrect')
        else
            fileName = varargin{1};
        end
        
    otherwise
        error('Please specify a valid analysis')
end
    
%% READ INPUT DATA
%  Input arguments passed from the system prompt will be received as strings
%  Thus, converting strings to double if required 

% X - expression data n-by-p
exist X 'var'
if ans
    if ~isa(X,'double')
        Xfile = X;     
        X = load(Xfile, '-ascii');
    end
end
 
%E - interaction matrix (network edges)
exist E 'var'
if ans
    if ~isa(E,'double')
        Efile = E;    
        E = load(Efile, '-ascii');
    end
end


%class - Y
exist Y 'var'
if ans
    if ~isa(Y,'double')
        Yfile = Y;
        Y = load(Yfile, '-ascii');
    end
end

whos lambda
% %lambdas for cross validation
% exist lambdas 'var'
% if ans
%     if ~isa(lambdas, 'double')
%         lambdas = str2num(lambdas);   
%     end
% end


% %folds for cross validation
% exist folds 'var'
% if ans
%     if ~isa(folds, 'double')
%         folds = str2num(folds);
%     end
% end

% Feature annotation
exist features 'var'
if ans
    if ~isa(features,'cell')
        featureAnnotfile = features;           
        features = table2cell(readtable(featureAnnotfile, 'ReadVariableNames', false));
    end
end

% Sample annotation
exist samples 'var'
if ans
    if ~isa(samples,'cell')
        samplesFile = samples;           
        samples = table2cell(readtable(samplesFile, 'ReadVariableNames', false));
    end
end




%% ADD PATH TO LIBRARIES
if ~isdeployed
    addpath('Libraries/')
    addpath('Libraries/fastGapFill/')
    addpath('Libraries/gaimc/')
    addpath('Libraries/GenePEN/')
    addpath('Libraries/TFOCS-1.3.1/')
end


%% Analysis

switch analysisType
    case 'hierarchicalClustering'
        % test
        samples={'sample1' 'sample2' 'sample3' 'sample4'};
        features={'f1' 'f2' 'f3' 'f4' 'f5' 'f6' 'f7'};
        hierarchicalClustering(X(1:4,1:7), samples, features)
    
    case 'crossValidation'                
        %crossValidation for different lambdas
        fprintf('Performing cross validation... \n')        
        [~, ~,stats, yTest, yTestPred] = crossValidation(X, E, Y, lambdas, folds, numIter);
        
        if ~strcmp(saveResults,'false')
            if strcmp(saveResults, 'true')
                outputDir = 'output_MultiPEN/crossValidation/';
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
        
          
    case 'featureSelection'
        %Feature selection for a specific lambda
        fprintf('Performing feature selection... \n')
        %FS is a table with columns: [name, weight, ranking]
        [FS, vts] = featureSelection(X, E, Y, features, lambda, numIter);
        
        %compute fold change
        [FS, higherControl, higherCases] = foldChange(FS, X, Y, samples);
        
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
            save([fileName '.mat'], 'FS')
            
            %intercept learnt from feature selection
            fileName = [outputDir 'MultiPEN-vts_lambda' num2str(lambda)];            
            dlmwrite([fileName '.txt'], vts, 'delimiter', '\t');    
            
            %% write separate tables for higherControl an higherCases
            fileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) '_higher-in-control.txt'];
            fprintf('Writing feature selection to file: \n\t%s\n',fileName)
            writetable(higherControl, fileName, 'delimiter', '\t');
            
            fileName = [outputDir 'MultiPEN-Rankings_lambda' num2str(lambda) '_higher-in-cases.txt'];
            fprintf('Writing feature selection to file: \n\t%s\n',fileName)
            writetable(higherCases, fileName, 'delimiter', '\t');
            %
        end
                       
        MP = FS;
        
        
    case 'enrichmentGO'
        
        % output directory 
        if strcmp(saveResults, 'true')
            outputDir = 'output_MultiPEN/enrichment-GO/';
        else
            outputDir = saveResults;
        end

        %check if output directory exists
        if exist(outputDir, 'dir') ~= 7
            mkdir(outputDir)
        end
        
        % Build string to call the R string enrichmentGO.R  (using Rscript)
        % syntaxis:
        % path_to_Rscript script_to_run file_name output_directory
        callToRscript = '/Library/Frameworks/R.framework/Resources/Rscript scripts/pathwayAnalysis/enrichmentGO.R';
        callToRscript = [callToRscript ' ' fileName ' ' outputDir];
        system(callToRscript)
        
        % Load the table with results
        MP = readtable([outputDir 'enrichment-GO.txt'], 'Delimiter', '\t');
        

end
function hierarchicalClustering(expression, samples, features, saveFigure, varargin)

% Plots clustergram for the expression data
% Inputs:
%   expression  Expression data, n-by-p, n samples, p features 
%   samples     Sample Annotation, n-by-1 cell vector
%   features    Feature Annotation, p-by-1 cell vector
%   saveFigure  Save figure?: 'true' 'false' outputDirectory
%   varargin 1  Threshold to filter expression data (exclusive)
%            2  Title plot

if numel(varargin)>=1
    %threshold
    threshold = varargin{1};
    ind = min(expression)>threshold;
    expression = expression(:,ind);
    features = features(ind,:);
    
    if numel(varargin)>1
        %Title Plot
        plotTitle = varargin(2); 
    else
        plotTitle = 'Hierarchical Clustering'; %default plot title
    end
else
    plotTitle = 'Hierarchical Clustering'; %default plot title
end


% Plot hierarchical clustering
% Standardising along the rows (features) before hierarchical clustering
% as expression transposed is p-by-n 
hc = clustergram(expression',...
    'ColumnLabels',samples,...    
    'Standardize',2,'Colormap','redbluecmap', ...
    'OptimalLeafOrder', true);
set(hc,'RowLabels',features)
addTitle(hc, plotTitle)


if ~strcmp(saveFigure,'false')
    close all hidden
    plot(hc)
    
    if strcmp(saveFigure, 'true')
        outputDir = 'output_MultiPEN/stats/';
    else
        outputDir = saveResults;
    end

    %check if output directory exists
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir)
    end

    % Save figure
    fileName = [outputDir 'hierarchical_clustering.png'];
    saveas(gcf, fileName)
    fprintf('\tFigure save to file: \n\t%s\n', ...
        fileName)
                        
end
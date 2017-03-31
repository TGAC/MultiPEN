function pcaExpressionData(X, groupLabels, titleGraph, annot, outputDir)
% Performs PCA
% Plots first three and two principal components
%
% Inputs:
% X is the n-by-p data matrix for n samples, p features
% groupLabels
% Example of zone's Labels
% labelZones=[{'PZ'}, {'PZ'}, {'PZ'}, {'PZ'}, {'PZ'}, {'PZ'}, {'PZ'}, ...
%    {'TZ'}, {'TZ'}, {'TZ'}, {'TZ'}, {'TZ'}, {'TZ'}, {'TZ'}];
% titleGraph
% annot     Annotation for features (genes/metabolites)
% outputDir  Location to save the figures

[coeff, scores, pcvars] = pca(X);

%pick the principal components
% x,y,z coordinates
x = zscore(scores(:,1));
y = zscore(scores(:,2));
z = zscore(scores(:,3));

figure()
c = unique(groupLabels);  % c contains the unique zones
switch numel(c)
    case 2
        gscatter3(x,y,z,groupLabels, {'b', 'm'}, {'.','.'}, 20);
    case 4 
        gscatter3(x,y,z,groupLabels, {'b', 'k', 'm', 'g'}, {'.','.','.','.'}, 20);
    case 10
        gscatter3(x,y,z,groupLabels, {'b', 'b', 'k', 'k', 'm', 'm', 'g', 'g', 'r', 'r'}, ...
            {'*','s','*','s','*','s','*','s','*','s'});
end
title(titleGraph)

xlabel(['PC-1(' num2str(round(pcvars(1)/sum(pcvars)*100)) '% variance)'])
ylabel(['PC-2(' num2str(round(pcvars(2)/sum(pcvars)*100)) '% variance)'])
zlabel(['PC-3(' num2str(round(pcvars(3)/sum(pcvars)*100)) '% variance)'])

if ~strcmp(outputDir,'false')
    % Check if output directory exists
    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir)
    end

    % Save figure
    fileName = [outputDir 'pca_three_components.png'];
    saveas(gcf, fileName)
    fprintf('\tFigure save to file: \n\t%s\n', ...
        fileName)
end

%% 
% Plot first two components

% plot the first two principal components
figure()
switch numel(c)
    case 2
        gscatter(x,y,groupLabels', 'bm', '*s');
    case 4
        gscatter(x,y, groupLabels', 'bkmr', '*sod')
    case 10
        gscatter(x,y,groupLabels', 'bbkkmmggrr', '*s*s*s*s*s');
end

xlabel(['PC-1(' num2str(round(pcvars(1)/sum(pcvars)*100)) '% variance)'])
ylabel(['PC-2(' num2str(round(pcvars(2)/sum(pcvars)*100)) '% variance)'])
title(titleGraph)
        
        
%for only the first two groups
% indx = find(strcmp(groupLabels, c(1)));
% plot(scores(indx,1),scores(indx,2),'+m')
% indx = find(strcmp(groupLabels, c(2)));
% hold on, plot(scores(indx,1),scores(indx,2),'+b'), hold off
% legend(c)
% xlabel(['PC-1(' num2str(round(pcvars(1)/sum(pcvars)*100)) '% variance)'])
% ylabel(['PC-2(' num2str(round(pcvars(2)/sum(pcvars)*100)) '% variance)'])
% title(titleGraph)

% Save figure
if ~strcmp(outputDir,'false')
    fileName = [outputDir 'pca_two_components.png'];
    saveas(gcf, fileName)
    fprintf('\tFigure save to file: \n\t%s\n', ...
        fileName)
end


%% 
%if annotation for each feature (gene/metabolite) is provided
% biplot of the pca coefficients 
if nargin > 3
    %biplot
    figure()
    biplot(coeff(:,1:3),'scores',scores(:,1:3),'varlabels',annot);
    title(titleGraph)
    
    % Save figure
    if ~strcmp(outputDir, 'false')
        fileName = [outputDir 'pca_biplot.png'];
        saveas(gcf, fileName)
        fprintf('\tFigure save to file: \n\t%s\n', fileName)
    end
end
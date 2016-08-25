function [Xprime, E, Y, features, sampleAnnot] = combinedNetwork(X, XAnnotation, XSamples, gClass)
%% Data Integration: transcriptomics and metabolomics - logistic regression
% X is the combined matrix n-by-p, where p = g+m (g is the number of genes
%    and m is the number of metabolites)
% XAnnotation is a p-by-1 cell vector with the geneSymbols and Chemical IDs


%% Load the network for the combined data
% loads table network [
load('../EscapePilot/MultiPEN/inputMultiPEN/combined/PC-escape_edgeList.mat')
network.weight = ones(numel(network.nodeA),1);  %since weigths are missing assign 1


%% list of unique nodes in the network

%remove edges which node(s) are empty
indx = find(strcmp(network.nodeA, ''));
if ~isempty(indx)
    network(indx,:) = [];
end

indx = find(strcmp(network.nodeB, ''));
if ~isempty(indx)
    network(indx,:) = [];
end

% create table for all unique nodes
nodes = table();
IDs = [network.nodeA; network.nodeB];
IDs = unique(IDs);
nodes.IDs = unique(IDs);
clear IDs


%% Build X prime for ONLY genes/metabolites that are in the network
%  and XAnnotation prime 

nodesUnique = {};
Xprime = [];
% check what elements in XAnnotation are in the 'nodes' list
for ii = 1:numel(XAnnotation)
%for ii = 1:10
    current = XAnnotation{ii};
    index = find(strcmp(current, nodes.IDs));
    if ~isempty(index)
        nodesUnique{end+1,1} = current;
        Xprime(:,end+1) = X(:,ii);
        %current
        %ii
    %else
        %fprintf('It is not on the nodes list')
    end
end
clear ii current index
% List of filtered data is now on Xprime and nodesUnique

%% Build E (indexed network prime) with indexes according to nodesUnique order
nodesUnique = cell2table(nodesUnique);
E = [];
for ii = 1: numel(network(:,1))
    current = [network.nodeA(ii) network.nodeB(ii)];
    index1 = find(strcmp(current{1}, table2array(nodesUnique)));
    index2 = find(strcmp(current{2}, table2array(nodesUnique)));
    if ~isempty(index1) && ~isempty(index2)
        E(end+1,:) = [index1 index2 1]; 
    end
end
clear ii current index


%% Write to file
% E - indexed edge list
outputDir = 'scripts/2016-08-22_Test-for_GSE46300_MTBLS174/MultiPEN_input-data/';
dlmwrite([outputDir 'E.txt'], E, 'delimiter', '\t')
% network - edge list with gene IDs and chEBI IDs
writetable(network, [outputDir 'edges.txt'], 'delimiter', '\t')
% combined data X - gene expression and metabolite levels
dlmwrite([outputDir 'X.txt'], Xprime, 'delimiter', '\t')
% combined data annotation
nodesUnique.Properties.VariableNames = {'XAnnotation'};
writetable(nodesUnique, [outputDir 'XAnnotation.txt'], 'Delimiter', '\t')
features = table2array(nodesUnique);
% XSamples - annotation for samples
XSamples = cell2table(XSamples);
XSamples.Properties.VariableNames = {'gSamples' 'mSamples'};
writetable(XSamples, [outputDir 'XSampleAnnotation.txt'], 'delimiter', '\t')
sampleAnnot = XSamples.gSamples;
% Y - low is used as control (zero value) and high as case (1)
% use column 1 which is gClass
Y = strcmp(gClass, 'high');
dlmwrite([outputDir 'Y.txt'], Y, 'delimiter', '\t')
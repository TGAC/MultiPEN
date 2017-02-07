function [E, networkT] = subnetwork4ExpressionData(network, XAnnotation)
% Function to create the subnetwork with only nodes that correspond to the 
% features for the expression data.

% network   The full network. It's a table with three columns: source, target and weight
% XAnnotation The p-by-1 cell vector with the feature's names 
%       For genes this is the geneSymbols, for metabolites is the Chemical IDs
%       preceded by the string 'CHEBI:', e.g. 'CHEBI:10049'

%% list of unique nodes in the network

% Remove edges which node(s) are empty
indx = find(strcmp(network.source, ''));
if ~isempty(indx)
    network(indx,:) = [];
end

indx = find(strcmp(network.target, ''));
if ~isempty(indx)
    network(indx,:) = [];
end

% create table for all unique nodes in network
nodes = table();
IDs = [network.source; network.target];
IDs = unique(IDs);
nodes.IDs = unique(IDs);
clear IDs


%% Expression Data 
%  does not contain duplicates
nn = numel(XAnnotation);
mm = numel(unique(XAnnotation));
if nn ~= mm
    error('There are duplicated feature names')
end


%% Build E 
% Nodes are represented by the index according to XAnnotation 
E = [];
networkT = table();
for ii = 1: numel(network(:,1))
    currentS = network.source{ii};
    currentT = network.target{ii};
    w = network.weight(ii);
    index1 = find(strcmp(currentS, XAnnotation));
    index2 = find(strcmp(currentT, XAnnotation));
    if ~isempty(index1) && ~isempty(index2)
        E(end+1,:) = [index1 index2 w]; 
        networkT(end+1,:) = network(ii,:);
    end
end

%% Display information
fprintf('Original Network has %i interactions \n', numel(network.source))
fprintf('of which %i include interactions for the expression data \n', numel(E(:,1)))
fprintf('\n')
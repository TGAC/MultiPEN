function [E, Xprime, XAnnotationPrime] = subnetwork4ExpressionData(network, X, XAnnotation)
% Function to create the subnetwork with only nodes that correspond to the 
% features for the expression data.

% network   The full network. It's a table with three columns: nodeA, nodeB and weight
% X     The expression matrix n-by-p, where p is the number of features
%       (genes and/or metabolites)
% XAnnotation The p-by-1 cell vector with the feature's names 
%       For genes this is the geneSymbols, for metabolites is the Chemical IDs

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


%% Build Xprime and XAnnotationPrime for ONLY the features 
%  (genes/metabolites) that are in the network 
XAnnotationPrime = {};
Xprime = [];
% check what elements in XAnnotation are in the 'nodes' list
for ii = 1:numel(XAnnotation)
    current = XAnnotation{ii};
    index = find(strcmp(current, nodes.IDs));
    %check that the current feature is not already on the list
    if isempty(find(strcmp(current, XAnnotationPrime)))
        if ~isempty(index)
            XAnnotationPrime{end+1,1} = current;
            Xprime(:,end+1) = X(:,ii);            
        end
    else
        error('There are duplicated feature names')
    end
end
clear ii current index
% List of filtered data is now on Xprime and XAnnotationPrime

%% Build E (indexed network prime) with indexes according to nodesUnique order
E = [];
for ii = 1: numel(network(:,1))
    current = [network.nodeA(ii) network.nodeB(ii)];
    w = network.weight(ii);
    index1 = find(strcmp(current(1), XAnnotationPrime));
    index2 = find(strcmp(current{2}, XAnnotationPrime));
    if ~isempty(index1) && ~isempty(index2)
        E(end+1,:) = [index1 index2 w]; 
    end
end
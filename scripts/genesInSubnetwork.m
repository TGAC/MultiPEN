function [XinNetwork, XAnnotationT] = genesInSubnetwork(X, XAnnotation, networkT)
% Get only expression of genes which are included in the network

%% First, compute subnetwork for expression data
%N = readtable(edges, 'delimiter', '\t');
%[E, networkT] = subnetwork4ExpressionData(N, X.name);

%% Extract unique genes in data
genes = XAnnotation;
fprintf('Number of unique genes: %i\n', numel(genes))


%% Extract unique nodes in Network
nodes = {networkT.source{:}, networkT.target{:}};
nodes = nodes';
nodes = unique(nodes);
fprintf('Number of unique nodes: %i\n', numel(nodes))


%% Build X 
% Nodes are represented by the index according to XAnnotation 
XinNetwork = [];
XAnnotationT = {};
for ii = 1: numel(genes)
    current = genes{ii};
    index = find(strcmp(current, nodes));    
    if ~isempty(index)
        XinNetwork(:,end+1) = X(:,ii);
        XAnnotationT{end+1,1} = current;
    end
end

%% Display information
fprintf('Original gene expression has %i genes \n', numel(genes))
fprintf('of which %i are nodes in the network \n', numel(XinNetwork(1,:)))
fprintf('\n')
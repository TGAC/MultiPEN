function E = subnetwork4ExpressionData(network, XAnnotation)
% Function to create the subnetwork with only nodes that correspond to the 
% features for the expression data.

% network   The full network. It's a table with three columns: source, source and weight
% X     The expression matrix n-by-p, where p is the number of features
%       (genes and/or metabolites)
% XAnnotation The p-by-1 cell vector with the feature's names 
%       For genes this is the geneSymbols, for metabolites is the Chemical IDs

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

% % %% Build Xprime and XAnnotationPrime for ONLY the features 
% % %  (genes/metabolites) that are in the network 
% % XAnnotationPrime = {};
% % Xprime = [];
% % % check what elements in XAnnotation are in the 'nodes' list
% % for ii = 1:numel(XAnnotation)
% %     current = XAnnotation{ii};
% %     index = find(strcmp(current, nodes.IDs));
% %     %check that the current feature is not already on the list
% %     if isempty(find(strcmp(current, XAnnotationPrime)))
% %         if ~isempty(index)
% %             XAnnotationPrime{end+1,1} = current;
% %             Xprime(:,end+1) = X(:,ii);            
% %         end
% %     else
% %         error('There are duplicated feature names')
% %     end
% % end
% % clear ii current index
% % % List of filtered data is now on Xprime and XAnnotationPrime

%% Build E 
% Nodes are represented by the index according to XAnnotation 
E = [];
for ii = 1: numel(network(:,1))
    currentS = network.source{ii};
    currentT = network.target{ii};
    w = network.weight(ii);
    index1 = find(strcmp(currentS, XAnnotation));
    index2 = find(strcmp(currentT, XAnnotation));
    if ~isempty(index1) && ~isempty(index2)
        E(end+1,:) = [index1 index2 w]; 
    end
end
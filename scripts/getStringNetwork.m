function gene2geneNetwork = getStringNetwork(geneName, index, stringID, network)
% Builds gene2geneNetwork with five columns:
%   from            ensembl protein name (used by StringDB)
%   to              ensembl protein name (used by StringDB)
%   score           Interaction score
%   fromIndex       Node's index in expression matrix 
%   toIndex         Node's index in expression matrix


% Inputs:
%   geneName        String Vector with gene names, p-by-1 cell array
%   index           Gene indexes in expression matrix, p-by-1 vector
%   stringID        StringDB ID for each gene, p-by-1 cell array
%   network         edge table for network: 
%                   [nodeA_EnsemblID, nodeB_ensembleID, score]
 
%count the number of genes with duplicated short names
genesWithDuplicatedShortName = 0;
 
%map all the edges of the network table
[nEdges,~] = size(network);
from = {};
to = from;
score = [];
fromIndx = [];
toIndx = [];
display(['Mapping ' num2str(nEdges) ' nodes'])
for ii = 1:nEdges
    if mod(ii,10000)==0, display(['    node ' num2str(ii) '...']), end    
    node1 = network.from(ii);
    node2 = network.to(ii);
    s = network.score(ii);
    %node's short names
    %subsetNode1 = mapped(strcmp(mapped.stringId, node1),:);
    position = find(strcmp(stringID, node1));
    subsetNode1 = table(geneName(position), stringID(position), ...
        'VariableNames', {'geneName' 'stringID'}); 
    %subsetNode2 = mapped(strcmp(mapped.stringId, node2),:);
    position = find(strcmp(stringID, node2));
    subsetNode2 = table(geneName(position), stringID(position), ...
        'VariableNames', {'geneName' 'stringID'});
    [duplicatesNode1,~] = size(subsetNode1);   
    [duplicatesNode2,~] = size(subsetNode2);      
    for jj = 1:duplicatesNode1        
        for kk = 1:duplicatesNode2
            %fromShortName = subsetNode1.shortName(jj);
            %toShortName = subsetNode2.shortName(kk);
            fromShortName = subsetNode1.geneName(jj);
            toShortName = subsetNode2.geneName(kk);
            %node's indexes
            %what's the index of each node (from and to)
            %from table geneIndex            
            %indexes to the tables: GenesNonZero and ExpressionNonZero
            %indxNode1 = geneIndex(strcmp(geneIndex.shortName, fromShortName),:);
            %indxNode2 = geneIndex(strcmp(geneIndex.shortName, toShortName),:);
            indxNode1 = index(strcmp(geneName, fromShortName));
            indxNode2 = index(strcmp(geneName, toShortName),:);
            [numberIndxNode1,~] = size(indxNode1);
            [numberIndxNode2,~] = size(indxNode2);
            if numberIndxNode1>1, genesWithDuplicatedShortName = genesWithDuplicatedShortName + 1; end
            if numberIndxNode2>1, genesWithDuplicatedShortName = genesWithDuplicatedShortName + 1; end
            for dd1 = 1: numberIndxNode1
                for dd2 = 1: numberIndxNode2                    
                    from(end+1) = fromShortName;
                    to(end+1) = toShortName;
                    score(end+1) = s;
                    %fromIndx(end+1) = indxNode1.index(dd1); 
                    %toIndx(end+1) = indxNode2.index(dd2);
                    fromIndx(end+1) = indxNode1(dd1); 
                    toIndx(end+1) = indxNode2(dd2);
                    % to validate this procedure I print:
                    % the current edge to be mapped
                    % the mapped shortName and indexes for node1 and node2
                    if dd1>1 || dd2>1
                        display(network(ii,:))
                        display(subsetNode1)
                        display(subsetNode2)
                        display(indxNode1)
                        display(indxNode2)
                        display([from(end-1:end) to(end-1:end) num2str(score(end-1:end)) ...
                            num2str(fromIndx(end-1:end)) num2str(toIndx(end-1:end))])
                        input('')
                    end
                end
            end              
        end
    end
end
 
fprintf('%i genes were found in multiple locations \n', genesWithDuplicatedShortName)
from = from';
to = to';
score = score';
fromIndx = fromIndx';
toIndx = toIndx';
gene2geneNetwork = table(from, to, score, fromIndx, toIndx, ...
    'VariableNames', {'from', 'to', 'score', 'fromIndx', 'toIndx'});
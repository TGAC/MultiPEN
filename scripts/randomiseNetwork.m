%function results = randomiseNetwork(X, E, elementToChange, typeOfChange, percentage)
function randomiseNetwork(X, E, Y, geneNames, lambda, outputDir, ...
    geneIndex, elementToChange, typeOfChange, percentage, numRandomNetworks, numIter)
%Randomise the network according to the type of change specified
% It generates a number of randomly modified networks (numRandomNetwork) 
% 

%   gene2geneNetwork
%example with network from escape:
%load('/Users/troncosp/Documents/Projects/EscapePilot/inputForGenePEN/gene2geneNetwork.mat')

rng('default')  % for reproducibility

%get general features
[numberEdges,~] = size(E);
numberNodes = numel(unique([E(:,1); E(:,2)]));
fprintf('ORIGINAL NETWORK FEATURES \n')
fprintf('Number of Edges: %i\n', numberEdges)
fprintf('Number of Nodes: %i\n', numberNodes)
%TO DO: degree distribution

switch elementToChange
    case 'edge'
        totalElements = numberEdges;
    case 'node'
        totalElements = numberNodes;
end

changedElements = [];
selectedGenes = [];
outputRandomNetwork = [outputDir 'randomEdges_originalNetwork/'];
geneSelectionE = GenePEN_test(X, E, Y, geneNames, lambda, outputRandomNetwork, geneIndex, numIter);

for ii = 1 : numRandomNetworks
    modifiedE = E;
    nChanges = ceil(totalElements*percentage);
    changes = randperm(totalElements,nChanges);

    switch typeOfChange
        case 'swap'         
            newToNodes = changes(end:-1:1);        
            modifiedE(changes,2) = E(newToNodes,2);
        case 'delete'
            fprintf('In development...')    
            modifiedE(changes,:) = [];
    end

    %where to save the random network
    outputRandomNetwork = [outputDir 'randomEdges_network' num2str(ii) '/'];
            
            
    %modifiedE is matrix [(from)index, (to)index, score]  
    changedElements(end+1,:) = changes;
%     save([outputRandomNetwork 'E.mat'], 'modifiedE')       
%            dlmwrite([outputRandomNetwork 'E.csv'], modifiedE)
    numIter = 300;
    geneSelection = GenePEN_test(X, modifiedE, Y, geneNames, lambda, outputRandomNetwork, geneIndex, numIter);
%             selectedGenes(:,end+1) = geneSelection.weight;
%             find(geneSelection.weight)'
end



for ii = 1 : numRandomNetworks
    fprintf('Modified Network, iteration %i\n', ii)
    fprintf('Changed Elements: \n')
    display(changedElements(ii,:))
%             fprintf('Selected Genes: \n')
%             display(find(selectedGenes(:,ii))')
end

display(unique(changedElements)')
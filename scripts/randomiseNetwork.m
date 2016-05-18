function randomNetwork = randomiseNetwork(E, elementToChange, change, percentage)
%Randomise the network according to the type of change specified
%   gene2geneNetwork

%example with network from escape:
%load('/Users/troncosp/Documents/Projects/EscapePilot/inputForGenePEN/gene2geneNetwork.mat')

%get general features
[numberEdges,~] = size(E);
numberNodes = numel(unique([E(:,1); E(:,2)]));
fprintf('ORIGINAL NETWORK FEATURES \n')
fprintf('Number of Edges: %i\n', numberEdges)
fprintf('Number of Nodes: %i\n', numberNodes)


switch elementToChange
    case 'edge'
        totalElements = numberEdges;
    case 'node'
        totalElements = numberNodes;
end

nChanges = ceil(totalElements*percentage);

switch change
    case 'random'        
        p = randperm(totalElements,nChanges);
        newToNodes = p(end:-1:1);
        from = E(:,1);        
        score = E(:,3);                        
        to =  E(:,2);
        to(p) = E(newToNodes,2);
        randomNetwork = [from to score];
    case 'leaf'
    case 'hub'
end




% function randomNetwork = randomizeNetwork(network, elementToChange, change, percentage)
% %Randomise the network according to the type of change specified
% %   gene2geneNetwork
% 
% %example with network from escape:
% %load('/Users/troncosp/Documents/Projects/EscapePilot/inputForGenePEN/gene2geneNetwork.mat')
% 
% %get general features
% [numberEdges,~] = size(network);
% numberNodes = numel(unique([network.fromIndx; network.toIndx]));
% fprintf('ORIGINAL NETWORK FEATURES \n')
% fprintf('Number of Edges: %i\n', numberEdges)
% fprintf('Number of Nodes: %i\n', numberNodes)
% 
% 
% switch elementToChange
%     case 'edge'
%         totalElements = numberEdges;
%     case 'node'
%         totalElements = numberNodes;
% end
% 
% nChanges = ceil(totalElements*percentage);
% 
% switch change
%     case 'random'        
%         p = randperm(totalElements,nChanges);
%         newToNodes = p(end:-1:1)
%         from = network.from;
%         fromIndex = network.fromIndx;        
%         score = network.score;        
%         to = network.to;
%         to(p) = network.to(newToNodes);
%         toIndex = network.toIndx;
%         toIndex(p) = network.toIndx(newToNodes);
%         randomNetwork = table(from,to, score, fromIndex, toIndex, ...
%             'VariableNames', {'from' 'to' 'score' 'fromIndx' 'toIndx'});
%     case 'leaf'
%     case 'hub'
% end

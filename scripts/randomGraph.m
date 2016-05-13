%Generate a random graph for:
%       n  number of nodes 
%       p  probability
%       d  node degree

%It uses Erdos-Renyi Random Graph
%  http://www.mathworks.com/matlabcentral/fileexchange/4206-erdos-renyi-random-graph

addpath 'ErdosRenyi_RandomGraph/'


% example: G = erdosRenyi(40,0.01,4);
G = erdosRenyi(n,p,d);

plotGraphBasic(G,6,1);
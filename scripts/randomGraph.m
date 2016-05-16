%Generate a random graph for:
%       n  number of nodes 
%       p  probability
%       d  node degree

%Using "Erdos-Renyi Random Graph" resource
%  http://www.mathworks.com/matlabcentral/fileexchange/4206-erdos-renyi-random-graph

addpath 'ErdosRenyi_RandomGraph/'

% example: G = erdosRenyi(40,0.01,4);
G = erdosRenyi(n,p,d);

plotGraphBasic(G,6,1);


%Using "MIT Network Analysis" resource
%  http://strategic.mit.edu/downloads.php?page=matlab_networks

addpath 'MIT_NetworkAnalysis/'
A = random_graph(5,0.5,10);  %adjacency matrix A
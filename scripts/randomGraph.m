function A = randomGraph(algorithm, N, p, d, E)
%Generate a random graph for:
%       N  number of nodes 
%       p  probability
%       d  node degree
%       E fixed number of edges

switch algorithm
    case 'RenyiErdos'
        %Using "Renyi-Erdos Random Graph" resource
        %  http://www.mathworks.com/matlabcentral/fileexchange/4206-erdos-renyi-random-graph
        %  path to function:
        addpath '~/Documents/MATLAB/ErdosRenyi_RandomGraph/'

        % example: G = erdosRenyi(40,0.01,4);
        G = erdosRenyi(N,p,d);
        A = G.Adj;

        plotGraphBasic(G,6,1);
    case 'MIT'
        %Using "MIT Network Analysis" resource
        %  http://strategic.mit.edu/downloads.php?page=matlab_networks

        addpath '~/Documents/MATLAB/MIT_NetworkAnalysis/'
        % random_graph(N, p, E, distribution, degrees)
         
        A = random_graph(N,p,E);  %adjacency matrix A       
end
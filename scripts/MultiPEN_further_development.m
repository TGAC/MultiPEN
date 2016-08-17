%MultiPEN: Further development

switch typeAnalysis
    case 'RandomiseNetwork'                
        %% UNDER DEVELOPMENT        
        outputDir = ['cuffnorm_output/string/GenePEN_results/' samplesForAnalysis '/' ];
        
        %load geneIndex
        load('cuffnorm_output/string/geneList_Index_1perRow.mat')        
                
        %elementToChange   can be 'node' or 'edge'
        %typeOfChange      can be 'swap', 'delete'
        %percentage of elements to change 
        %test with just first 100 genes and 
        %corresponding network (variable E100Nodes):
        load('cuffnorm_output/string/GenePEN_results/test2/E100Genes.mat')
        X100 = X(:,1:100);
        featureNames = featureNames(1:100,:);  
        elementToChange = 'edge';
        typeOfChange = 'swap';
        perc = .40;   
        numRandomNetworks = 5;
            
        %randomiseNetwork will run geneSelection which needs the following
        %inputs:
        %geneSelection(X, modifiedE, Y, geneNames, lambdas, outputRandomNetwork, geneIndex);
        
        randomiseNetwork(X100, E100Nodes, Y, featureNames, lambda, outputDir, ...
            geneIndex, elementToChange, typeOfChange, perc, numRandomNetworks);
        
    case 'ErdosRenyi'
        %% UNDER DEVELOPMENT 
        outputDir = ['cuffnorm_output/string/GenePEN_results/' samplesForAnalysis '/RenyiErdos/' ]; 
        %load geneIndex
        load('cuffnorm_output/string/geneList_Index_1perRow.mat')
        %test with just first 100 genes and 
        %corresponding network (variable E100Nodes):
        %load network
        %load('cuffnorm_output/string/GenePEN_results/test2/E100Genes.mat')
        X100 = X(:,1:100);
        geneNames100 = featureNames(1:100,:);
        numRandomNetworks = 5;       
        algorithm = 'MIT';
        %get the corresponding network for those 100 genes
        E100Nodes = [];
        for in100 = 1 : numel(E(:,1))
            indxs1 = find(geneNames100.index == E(in100,1));
            indxs2 = find(geneNames100.index == E(in100,2));                        
            if ~isempty(indxs1) && ~isempty(indxs2)
                E100Nodes(end+1,:) = [indxs1 indxs2 E(in100,3)]; 
            end
        end
        
        %run RenyiErdos test
        %   RenyiErdos test will run geneSelection which needs the following
        %inputs:
        %geneSelection(X, modifiedE, Y, geneNames, lambdas, outputRandomNetwork, geneIndex);
        
        ErdosRenyiNetwork(X100, E100Nodes, Y, geneNames100, lambda, outputDir, ...
            geneIndex, numRandomNetworks, algorithm);
end
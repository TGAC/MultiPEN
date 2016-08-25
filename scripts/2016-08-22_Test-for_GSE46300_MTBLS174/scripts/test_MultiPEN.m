%% LOAD GENE EXPRESSION DATA
[geneExpression, illuminaProbes, samples, sampleGEO] = loadData();

% Get the class per sample from the samples 
% 

%% save illuminaProbes for mapping to ensembl ID and gene symbols in R
illuminaProbes = cell2table(illuminaProbes);
illuminaProbes.Properties.VariableNames = {'ProbeID'};
writetable(illuminaProbes, 'MultiPEN_input-data/illuminaProbes.txt')
illuminaProbes = table2cell(illuminaProbes);

%% gene annotation file which is a table with columns:
% entrezID, probeID, gSymbols
fileName = '2016-08-22_Data-Analysis/MultiPEN_input-data/annot_entrezID_probeID_gSymbols.txt';
annot_gene = readtable(fileName, 'delimiter', '\t');

%order annotation according to illuminaProbes order
n = numel(illuminaProbes);
entrezID = cell(n,1);
gSymbol = cell(n,1);
for ii=1:n
    p = illuminaProbes{ii};
    indx = find(strcmp(annot_gene.probeID, p));
    if isempty(indx)
        entrezID{ii} = '';
        gSymbol{ii} = '';
    else
        entrezID{ii} = annot_gene.entrezID{indx};
        gSymbol{ii} = annot_gene.gSymbols{indx};        
    end
end



%% Load metabolites
%% Network for this dataset
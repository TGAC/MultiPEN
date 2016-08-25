%% LOAD GENE EXPRESSION DATA
[gExpression, illuminaProbes, gSamples, gSamplesGEO, gClass] = loadGeneExpression();

% Get the class per sample from the samples 
% class{1} contains a string of the type: 
%     grade of steatosis: low

for ii = 1:numel(gClass)
    gClass{ii} = strrep(gClass{ii}, 'grade of steatosis: ', '');
end



%% LOAD METABOLITES
[chemicalID, mSamples, mLevels, mIdentif] = loadMetabolomics();


%% SAVE illuminaProbes for mapping to ensembl ID and gene symbols in R
illuminaProbes = cell2table(illuminaProbes);
illuminaProbes.Properties.VariableNames = {'ProbeID'};
writetable(illuminaProbes, 'MultiPEN_input-data/illuminaProbes.txt')
illuminaProbes = table2cell(illuminaProbes);

%% LOAD MAPPING FOR ILLUMINA PROBES TO SYMBOLS AND ENSEMBL
%gene annotation file which is a table with columns:
% entrezID, probeID, gSymbols
fileName = 'scripts/2016-08-22_Test-for_GSE46300_MTBLS174/MultiPEN_input-data/annot_entrezID_probeID_gSymbols.txt';
annot_gene = readtable(fileName, 'delimiter', '\t');

%order annotation according to illuminaProbes order
n = numel(illuminaProbes);
entrezID = cell(n,1);
gSymbol = cell(n,1);
for ii=1:n
    p = illuminaProbes{ii};
    indxG = find(strcmp(annot_gene.probeID, p));
    if isempty(indxG)
        entrezID{ii} = '';
        gSymbol{ii} = '';
    else
        entrezID{ii} = annot_gene.entrezID{indxG};
        gSymbol{ii} = annot_gene.gSymbols{indxG};        
    end
end

clear fileName annot_gene n ii p



%% REMOVE DUPLICATED FOR gSymbols

%first remove all genes with no gene Symbol
indx = find(strcmp(gSymbol,''));
gSymbolT = gSymbol;
gSymbolT(indx) = [];  % remove genes with no gene Symbol
gExpressionT = gExpression;
gExpressionT(:,indx) = [];  % remove genes with no gene Symbol

uSymbols = unique(gSymbolT);  %unique gene symbols

nG = numel(uSymbols);  %number of unique gene symbols
gExpressionTU = zeros(numel(gExpressionT(:,1)), nG);
gAnnot = cell(nG,1);

for ii=1:nG
    g = uSymbols(ii);
    indx = find(strcmp(gSymbolT,g));
    
    if numel(indx) == 1
        gExpressionTU(:,ii) = gExpressionT(:,indx);
        gAnnot{ii} = g;
    else
        gExpressionTU(:,ii) = gExpressionT(:,indx(1));
        gAnnot{ii} = g;
        
        %delete the duplicates from gExpressionT and gSymbolT
        gExpressionT(:,indx) = [];
        gSymbolT(indx,:) = [];
    end
    %pause(1)
    %ii
    %indx
end

% Gene Expression and gene Symbols are in:
% gExpressionTU and gAnnot


%% Combined Data
% there are only eight samples with both 
% gene expression and metabolomics data
% H0004, H0007, H0008, H0009, H0011, H0018, H0021, H0022

% Gene Expression Data
% Remove H0012 sample as there is no metabolomics data
indxG = [];
for ii=1:numel(gSamples)
    if isempty(strfind(gSamples{ii}, '0012'))
        indxG(end+1,1) = ii;
    end
end
Xgenes = gExpressionTU(indxG,:);
gClassT = gClass(indxG);
%XAnnotation = gSamples(indxG);

% Metabolomics Data
% indxM = [];
% for ii = 1:numel(XAnnotation)
%     current = strsplit(XAnnotation{ii}, '_')
%     ID = current{2};
%     ID(1) = [];
%     find the indexes for the metabolomics samples
%     %find(isempty(strfind(mSamples, ID)))
% end
indxM = [17 17 16 16 15 15 13 13 18 18 14 14 7 7 8 8];




%% create the combined matrix
X = [Xgenes mLevels(indxM,:)];
XAnnotation = [gAnnot; chemicalID];
XSamples = [gSamples(indxG) mSamples(indxM)];



%% Network for this dataset
[Xnetwork, E, Y, features, sampleAnnot] = combinedNetwork(X, XAnnotation, XSamples, gClassT);
Y = double(Y);

%% Run MultiPEN - cross validation
pathMP = 'scripts/2016-08-22_Test-for_GSE46300_MTBLS174/MultiPEN_output-data/';
numIter = '100';
lambdas = '[0.00001 0.000055 0.0001 0.00055 0.001 0.0055 0.01 0.055 0.1 0.55 1]';
MP = MultiPEN('cross_validation', pathMP, Xnetwork, E, Y, lambdas, '1', numIter);

%%
lambda = '0.0001';
MP = MultiPEN('featureSelection', pathMP, Xnetwork, E2, Y, lambda, features, sampleAnnot, numIter);
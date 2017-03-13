function [FS, higherControl, higherCases] = foldChange(FS, X, Y, sampleAnnot)%, outputDir, samplesForAnalysis, expressionType)
%% Calculate fold change from 
%  to find up regulated and down-regulate features
%  FS is a p-by-3 table with columns: [name, weight, ranking]
%  X the expression data, a n-by-p table 

% ratio of the average between classes:  
% (mu(class1) / mu(class0)) - 1

nf = numel(FS.name);  %number of features
fchange = zeros(nf,1);
higherIn = repmat({''}, nf, 1);
higherControl = table();
higherCases = table();


for ii = 1:nf    
    indxControl = Y~=1;  %control - A
    control = X(indxControl, ii); %A
    
    indxCases = Y==1;  %case - B
    cases = X(indxCases, ii); %B
    
    % foldChange = (B/A) - 1
    r = (mean(cases) / mean(control)) - 1;
    fchange(ii) = r;
    % if foldChange is positive, B (case) is higher    
    if r > 0 
        higherIn(ii) = {'case'};
        % save list of higher in class 0 separately
        higherCases(end+1,:) = [FS(ii,:) higherIn(ii) array2table(X(:,ii)')];
    % if foldChange is negative, A (control) is higher    
    else 
        higherIn(ii) = {'control'};
        % save list of higher in class 1 separately
        higherControl(end+1,:) = [FS(ii,:) higherIn(ii) array2table(X(:,ii)')];
        
    end    
end

%add column names 
higherControl.Properties.VariableNames(4) = {'higherIn'};
higherControl.Properties.VariableNames(5:end) = sampleAnnot;
higherCases.Properties.VariableNames(4) = {'higherIn'};
higherCases.Properties.VariableNames(5:end) = sampleAnnot;


%% Add a column to feature selection 
FS.foldChange = fchange;
FS.higherIn = higherIn;
% FS(:,end+1:end+(numel(X(:,1)))) = array2table(X');
% FS.Properties.VariableNames(6:end) = sampleAnnot;
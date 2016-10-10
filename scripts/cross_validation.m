function [weights, vts, stats] = ...
    crossValidation(X, E, Y, all_lambdas, cvfold, numIter)

% Supervised sample classification and identification of
% network-grouped predictive features in gene and/or metabolite
% expression data using GenePEN method
% 
% This function runs cross validation on gene expression and/or 
% metabolite level data, and molecular interaction data for
% GenePEN 

% Inputs:
%   XDataFile       data matrix with rows = samples, columns = features
%                   (genes and or metabolites)
%                   assumes that data has been normalised/standardised
%   YDataFile       class for each sample, i.e., 0 for control, 1 for cases
%   EDataFile       Interaction matrix, i.e., the network's edge list. 
%                   Each row specifies an interaction in the format:
%                   [node node interaction-weight]
%   all_lambdas     select lambda parameters to be used
%                   (here either a single lambda parameter can be
%                   specified or multiple lambda parameters can be
%                   tested) e.g.:
%                   all_lambdas = logspace(-12,2,num_of_lambdas);  
%                   all_lambdas = logspace(-0.3,-0.01,num_of_lambdas);
%   cvfold          the number of cross-validation cycles
%   numIter         maximum number of itrations for optimisation
%       
% Outputs:
%   weights
%   vts  
%   stats           Number selected features, LCC (avg, std), AUC (avg,
%                   std), True Positive Rate


%% Show input parameters in screen
disp('lambdas for cross validation')
disp(all_lambdas)

disp('number of folds for cross validation: ')
disp(cvfold)

fprintf('Number of max iterations: \n\t%i\n\n', numIter)


%% Preparing Data
% store interaction matrix in sparse format
P = sparse(E(:,1),E(:,2),E(:,3),size(X,2),size(X,2),size(E,1));

% convert triangle matrix to symmetric matrix
A = P + P' - diag(diag(P));


%% Weights, vts and stats
% weights has m columns, where m = length(all_lambdas) * cvfold;
weights = zeros(numel(X(1,:)), length(all_lambdas)*cvfold );
vts = zeros(1, length(all_lambdas)*cvfold );
%statistics
AUC = zeros(cvfold,1);
LCC = AUC; 
selected = AUC;
avgAUC = zeros(length(all_lambdas),1);
stdAUC = avgAUC;
avgLCC = avgAUC;
stdLCC = avgAUC;
avgsel = avgAUC;

%% Compute cross-validation per lambda using GenePEN
for i = 1:length(all_lambdas)
    lambda = all_lambdas(i);
    fprintf('\n ######Lambda: %d\n', lambda)

    % use this line to ensure reproducibility    
    rng('default')
    
    c = cvpartition(Y,'kfold',cvfold);
    
    %Cross validation
    for j=1:cvfold
            %GenePEN
            Xtrain = X(c.training(j),:);
            Ytrain = Y(c.training(j),:);
            [wt, vt] = GenePEN( Xtrain, Ytrain, A, lambda, numIter ); 
            % index for the current column considering folds and lambdas
            indx = j + ((i-1)*cvfold) ;  
            weights(:,indx) = wt;   
            vts(1,j) = vt;       
            weights(abs(wt)<1e-8,indx) = 0;
            
            % Selected Features (indexes)
            S = find(abs(wt)>1e-8);

            % LCC's size
            if (size(S, 1) > 0)
                [~,p] = largest_component(A(S,S));
            else
                p = 0;
            end

            % Evaluate Performance
            Xtest = X(c.test(j),:);
            y = Y(c.test(j),:);
            y_pred = Xtest *  wt;

            if (max(y_pred) - min(y_pred) == 0)
              if(max(y_pred) < 0)
                ypred = zeros(1, size(y_pred,1));
              else
                ypred = ones(1, size(y_pred,1));
              end
           else
              ypred = (y_pred - min(y_pred))/(max(y_pred) - min(y_pred));
           end

           r = tiedrank(ypred);
           %AUC
           AUC(j) = (sum(r(y==1)) - sum(y==1) * (sum(y==1)+1)/2) /( sum(y<1) * sum(y==1)); 
           LCC(j) = sum(p); %size of the largest connected component
           selected(j) = numel(S); % selected features   
           
           disp(' ');
           disp('Current cycle:');
           disp(horzcat('Area under the curve (AUC): ',num2str(AUC(j))));
           disp(horzcat('Selected Features: ',num2str(selected(j))));
           disp(horzcat('Size of largest connected component (LCC): ',num2str(LCC(j))));
           disp(' ');       
           
           
    end

    % performance per lambda
    avgsel(i) = mean(selected);
    avgLCC(i) = mean(LCC);
    stdLCC(i) = std(LCC);
    avgAUC(i) = mean(AUC);
    stdAUC(i) = std(AUC);
    
    disp(' ');
    disp('-------------');
    disp(horzcat('Average AUC for selected lambda: ',num2str(avgAUC(i))));
    disp(' ');

end

% Cross Validation Performance
stats = table();
stats.lambda = all_lambdas';
stats.selected = avgsel;
stats.avgLCC = avgLCC;
stats.stdLCC = stdLCC;
stats.avgAUC = avgAUC;
stats.stdAUC = stdAUC;

display(stats)
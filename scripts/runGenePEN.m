function [weights, vt] = runGenePEN(X, E, Y, lambda, numIter)
% This function runs GenePEN [1] on expression data (gene expression 
% and/or metabolite level data) and a molecular interaction matrix.
% 
% [1] Vlassis N and Glaab E., Stat. Appl. Genet Mol Biol. 2015, 
%     doi: 10.1515/sagmb-2014-0045)

% Inputs:
%   X           expression matrix with rows = samples, columns = features
%               features are genes and/or metabolites
%               We assume that data has been normalised/standardised
%   E           Interaction matrix, i.e., the network's edge list. 
%               Each row specifies an interaction in the format:
%               [nodeA nodeB interaction-weight]
%   Y           sample class for each sample: 0 for control, 1 for cases
%   lambda      tradeoff parameter for the regularisation function
%   numIter     maximum number of iterations for optimisation
%       
% Outputs:
%   weights     weights calculated by GenePEN for each fearture
%   vts         the intercept
%   Ypred       predicted classes for samples in test set

%% Preparing Data
% store interaction matrix in sparse format
P = sparse(E(:,1),E(:,2),E(:,3),size(X,2),size(X,2),size(E,1));

% convert triangle matrix to symmetric matrix
A = P + P' - diag(diag(P));

[weights, vt] = GenePEN( X, Y, A, lambda, numIter ); 
weights(abs(weights)<1e-8) = 0;

% % (indexes) Selected Features 
% S = find(abs(weights)>1e-8);
% disp(horzcat('Selected Features: ',num2str(numel(S))));
% 
% %% Evaluate performance
% y_pred = X * weights;
% if (max(y_pred) - min(y_pred) == 0)
%     if(max(y_pred) < 0)
%         Ypred = zeros(1, size(y_pred,1));
%     else
%         Ypred = ones(1, size(y_pred,1));
%     end
% else
%     Ypred = (y_pred - min(y_pred))/(max(y_pred) - min(y_pred));
% end
% 
% % To evaluate prediction
% r = tiedrank(Ypred);
% %AUC
% auc = (sum(r(Y==1)) - sum(Y==1) * (sum(Y==1)+1)/2) /( sum(Y<1) * sum(Y==1)); 
% 
% % LCC's size
% if (size(S, 1) > 0)
%     [~,p] = largest_component(A(S,S));
% else
%     p = 0;
% end
% LCC = sum(p); %size of the largest connected component
% disp(horzcat('LCC: ',num2str(LCC)));
% disp(horzcat('auc: ',num2str(auc)));

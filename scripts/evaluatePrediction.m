function stats = evaluatePrediction(X, Y, E, weights, D)
% 
% D     decision threshold, typically 0.5 but for some application could be
%       different

%% Preparing Data
% store interaction matrix in sparse format
P = sparse(E(:,1),E(:,2),E(:,3),size(X,2),size(X,2),size(E,1));

% convert triangle matrix to symmetric matrix
A = P + P' - diag(diag(P));

%% Evaluate performance
y_pred = X * weights;
if (max(y_pred) - min(y_pred) == 0)
    if(max(y_pred) < 0)
        Ypred = zeros(1, size(y_pred,1));
    else
        Ypred = ones(1, size(y_pred,1));
    end
else
    Ypred = (y_pred - min(y_pred))/(max(y_pred) - min(y_pred));
end

% To evaluate prediction
r = tiedrank(Ypred);
%% AUC
auc = (sum(r(Y==1)) - sum(Y==1) * (sum(Y==1)+1)/2) /( sum(Y<1) * sum(Y==1)); 
disp(horzcat('auc: ',num2str(auc)));

%% LCC's size
S = find(abs(weights)>1e-8);
if (size(S, 1) > 0)
    [~,p] = largest_component(A(S,S));
else
    p = 0;
end
LCC = sum(p); %size of the largest connected component
disp(horzcat('LCC: ',num2str(LCC)));


%% Accuracy
% Accuracy
accuracy = sum(Y == (Ypred>D)) / numel(Y);
% True positive
indx = find(Y);
TP = sum( Y(indx) == (Ypred(indx) > D) ) / numel(indx);
% True negatives
indx = find(~Y);
TN = sum( Y(indx) == (Ypred(indx) > D) ) / numel(indx);
% False positive
indx = find(Y);
FP = sum( Y(indx) ~= (Ypred(indx) > D) ) / numel(indx);
% False negatives
indx = find(~Y);
FN = sum( Y(indx) ~= (Ypred(indx) > D) ) / numel(indx);

%% Stats
stats = table();
stats.LCC = LCC;
stats.auc = auc;
stats.accuracy = accuracy;
stats.TP = TP;
stats.TN = TN;
stats.FP = FP;
stats.FN = FN;
function vifTbl = tsglmm_compute_vif(mdl)

% Extract fixed-effects design matrix
X = mdl.designMatrix('Fixed');
names = mdl.CoefficientNames(:);

% Exclude intercept (first column)
interceptIdx = strcmp(names, '(Intercept)');
X(:,interceptIdx) = [];
names(interceptIdx) = [];

% Compute VIFs
P = size(X,2);
VIF = zeros(P,1);

for j = 1:P
    y = X(:,j);
    Xj = X(:, setdiff(1:P,j));
    % Add constant term
    Xj = [ones(size(Xj,1),1) Xj];
    % Regress y ~ others
    b = Xj\y;
    yhat = Xj*b;
    SSres = sum((y - yhat).^2);
    SStot = sum((y - mean(y)).^2);
    R2 = 1 - SSres/SStot;
    VIF(j) = 1/(1-R2);
end

% Return as table
vifTbl = table(names, VIF, 'VariableNames', {'Predictor','VIF'});
end

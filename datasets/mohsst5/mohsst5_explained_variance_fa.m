function v = mohsst5_explained_variance_fa(W, CovW, X, CovX)

% This is approximate, because this does not take into account the
% correlations between different time instances or spatial locations.

% Remove temporal mean
X = bsxfun(@minus, X, mean(X,2));

% Compute <XX>
XX = X*X';
if ndims(CovX) == 2 && all(size(CovX)==size(X))
  XX = XX + diag(sum(CovX,2));
else
  XX = XX + sum(CovX,3);
end

% Compute weighted <WW>
w = mohsst5_weights();
w = mohsst5_remove_land(w);
WW = W*diag(w)*W';
if ndims(CovW) == 2 && all(size(CovW)==size(W))
  WW = WW + diag(wsum(CovW,w,2));
else
  WW = WW + wsum(CovW,w,3);
end

% Effective number of samples
%N = sum(sum( w(:)*ones(1,size(X,2)) ));
N = sum(w) * size(X,2);

% Explained variance
v = diagprod(WW,XX) / N;
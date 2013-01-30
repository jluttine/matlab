% Expands singleton axes of the matrix.
%
% X = BROADCAST(X, D)
%
% D is the desired size. For non-singleton dimensions of X must hold
% SIZE(X,I) == D(I).

function X = broadcast(X, d)
X = bsxfun(@plus, zeros(d), X);

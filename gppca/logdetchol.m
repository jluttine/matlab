function logdetx = logdetchol(U)
% function logdetx = chollogdet(U)
%
% Returns log(det(X)) where X=U'*U and U is upper triangular matrix, that
% is, U is Cholesky decomposition of X.

logdetx = 2 * sum(log(diag(U)));
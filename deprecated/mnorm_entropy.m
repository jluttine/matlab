function H = mnorm_entropy(U)
warning('This function is deprecated')

% H = mnorm_entropy(U)
%
% U is the Cholesky factor of a covariance matrix C, that is, a lower or
% upper triangylar matrix such that C=U*U' or C=U'*U.

d = size(U,1);

H = d/2*log(2*pi) + 0.5*logdet_chol(U) + d/2;

% LOGDET - Computes the log-determinant of a square matrix X.
%
% Y = LOGDET(X)
%
% More precisely, the function takes the logarithm with respect to the
% absolute value of the determinant, i.e., Y = LOG(ABS(DET(X))). However,
% the function uses QR decomposition for stable and efficient performance.

% Last modified 2010-06-04
% Copyright (c) Jaakko Luttinen

function y = logdet(X)

[Q,R] = qr(X);
y = sum(log(abs(diag(R))));
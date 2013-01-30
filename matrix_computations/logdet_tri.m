% LOGDET_TRI - Computes the log-determinant of a triangular matrix.
%
% Y = LOGDET_TRI(U)
%
% U must be either lower or upper triangular square matrix.

% Last modified 2010-06-04
% Copyright (c) Jaakko Luttinen

function y = logdet_tri(U)

y = sum(log(diag(U)));
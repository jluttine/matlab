% LOGDET_CHOL  -  Computes the log-determinant of a matrix using its
% Cholesky factor.
%
% Y = LOGDET_CHOL(U)
%
% Computes the log-determinant of the matrix X whose Cholesky factor is
% the triangular matrix U. U can be lower or upper triangular.
%
% More generally, X could be a product of two instances of U and/or U'. For
% instance, X=U'*U' or X=U'*U.

% Last modified 2010-06-04
% Copyright (c) Jaakko Luttinen

function y = logdet_chol(U)

y = 2 * logdet_tri(U);
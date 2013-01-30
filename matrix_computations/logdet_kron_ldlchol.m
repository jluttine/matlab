% LOGDET_KRON_LDLCHOL - Computes the log-determinant of a matrix using its
% Cholesky factor from SuiteSparse.
%
% Y = LOGDET_LDLCHOL(LD1,LD2)
%
% Computes the log-determinant of the matrix X whose Cholesky factor is
% the triangular matrix U. U can be lower or upper triangular.
%
% More generally, X could be a product of two instances of U and/or U'. For
% instance, X=U'*U' or X=U'*U.

% Last modified 2010-11-09
% Copyright (c) Jaakko Luttinen

function y = logdet_ldlchol(LD1,LD2)

N1 = size(LD1,1);
N2 = size(LD2,1);

y = N2*logdet_tri(LD1) + N1*logdet_tri(LD2);
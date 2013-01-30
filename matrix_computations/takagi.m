% TAKAGI  -  Factorizes a positive semi-definite symmetric matrix.
%
% V = TAKAGI(A)
%
% The matrix A must be positive semi-definite and symmetric. The resulting
% matrix V is such that V*V'=A.  
%
% Note: This is not really Takagi decomposition but I couldn't think of a
% better name.. :(

% Last modified 2010-06-04
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function V = takagi(A)

[L,D] = ldl(A);
V = L * sqrt(D);

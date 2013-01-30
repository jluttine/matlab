% SQRT_COV  -  Factorizes a positive semi-definite symmetric matrix.
%
% V = SQRT_COV(A)
%
% The matrix A must be positive semi-definite and symmetric. The resulting
% matrix V is such that V*V'=A.  

% Last modified 2011-10-20
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function V = sqrt_cov(A)

[L,D] = ldl(A);
D(D<0) = 0;
V = L * sqrt(D);

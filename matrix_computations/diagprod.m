% DIAGPROD - Computes the diagonal elements of a matrix product.
%
% D = DIAGPROD(A,B)
%
% Evaluates diag(A*B') efficiently. A and B must have the same size.
%
% D = DIAGPROD(A)
%
% Evaluates diag(A*A') efficiently.

% Last modified 2010-06-24
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function d = diagprod(A,B)

if nargin < 2
  d = dot(A,A,2);
else
  d = dot(A,B,2);
end

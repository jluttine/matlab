function A = adddiag(A,v,k)
%ADDDIAG  Modify the diagonal(s) of a matrix.
% ADDDIAG(A,V)
% ADDDIAG(A,V,K)
%
% See also DIAG, FINDDIAG.

if nargin < 3
  k = 0;
end
i = finddiag(A,k);
A(i) = A(i) + v;

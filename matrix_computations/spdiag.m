% SPDIAG - Creates a sparse diagonal matrix from a given vector.
%
% D = SPDIAG(X)
%
% D is a sparse matrix that has the vector X on its main diagonal.
%
% D = SPDIAG(X,I)
%
% D is a sparse matrix that has the vector X on its I-th diagonal.

% Last modified 2010-06-18
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function D = spdiag(d,i)

if nargin < 2 || (isscalar(i) && i==0)
  m = numel(d);
  D = sparse(1:m, 1:m, d, m, m);
else
  m = numel(d) + abs(i);
  if i >= 0
    D = sparse(1:numel(d),(1:numel(d))+i,d,m,m);
  else
    D = sparse((1:numel(d))-i,(1:numel(d)),d,m,m);
  end
end


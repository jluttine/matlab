% LINSOLVE_CHOL - Solves C*X=Y using the Cholesky decomposition of C.
%
% A function call
% 
%   X = LINSOLVE_CHOL(R,Y)
%
% assumes C=R'R.
%
% A function call
%
%   X = LINSOLVE_CHOL(L,Y,'lower')
%
% assumes C=L*L'.
%
% R or L can be found using, e.g., the Cholesky decomposition CHOL.

% Last modified 2010-10-18
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Y = linsolve_lchol(L,Y,q)

% $$$ if nargin < 3
% $$$   x = ldlsolve(LD,x);
% $$$ else
% $$$   x(q(:),:) = ldlsolve(LD,x(q(:),:));
% $$$ end

if issparse(L)
  if nargin < 3
    Y = L' \ (L \ Y);
  else
    Y(q(:),:) = L' \ (L \ Y(q(:),:));
  end
else
  if nargin < 3
    opts.LT = true;
    Y = linsolve(L, full(Y), opts);
    opts.TRANSA = true;
    Y = linsolve(L, Y, opts);
  else
    error('Not yet implemented');
  end
end

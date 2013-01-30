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

function Y = linsolve_chol(R,Y,type)

if nargin >= 3 && strcmpi(type, 'lower')
  if issparse(R)
    Y = R' \ (R \ Y);
  else
    opts.LT = true;
    Y = linsolve(R, full(Y), opts);
    opts.TRANSA = true;
    Y = linsolve(R, Y, opts);
  end
else
  if issparse(R)
    Y = R \ (R' \ Y);
  else
    opts.UT = true;
    opts.TRANSA = true;
    Y = linsolve(R, full(Y), opts);
    opts.TRANSA = false;
    Y = linsolve(R, Y, opts);
  end
end

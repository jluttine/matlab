% LINSOLVE_LDLCHOL - Solves a matrix-vector equation when the matrix is a
%                    sparse symmetric positive definite matrix.
%
% Solves K*Y = X. Usage:
%
%   Y = LINSOLVE_LDLCHOL(LD,X)
%
% where LD = LDLCHOL(K).
%
% The function is just a simple wrapper in accordance with the naming
% conventions.
%
% See also LDLCHOL, LDLSOLVE.

% Last modified 2011-01-31
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = linsolve_ldlchol(LD,x,q)

%fprintf('Sparsity in linsolve_ldlchol: %f\n', sparsity(LD));

if nargin < 3
  x = ldlsolve(LD,x);
else
  x(q(:),:) = ldlsolve(LD,x(q(:),:));
end

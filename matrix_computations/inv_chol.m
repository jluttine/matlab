
function U = inv_chol(U, type)

if nargin < 2
  U = linsolve_chol(U,speye(size(U)));
else
  U = linsolve_chol(U,speye(size(U)),type);
end

% $$$ if nargin >= 2 && strcmpi(type, 'lower')
% $$$ else
% $$$   X = linsolve_triu(U,speye(size(U)));
% $$$   X = linsolve_triu(U,speye(size(U)));
% $$$ %  X = A*A';
% $$$ end


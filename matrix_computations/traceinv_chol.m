function y = traceinv_chol(U, type)

if nargin >= 2 && strcmpi(type, 'lower')
  A = linsolve_tril(U,speye(size(U)));
else
  A = linsolve_triu(U,speye(size(U)));
end

y = traceprod(A,A);
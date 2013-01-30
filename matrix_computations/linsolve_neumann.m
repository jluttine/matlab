% Approximates X = INV(A)*Y by a truncated series.
function x = linsolve_neumann(A,y,D,n)

if nargin < 3 && isnumeric(A)
  tol = 1e-6;
  D = (1+tol)*normest(A,tol); % this is probably veeeery bad...
end
if nargin < 4
  n = 100;
end

if isnumeric(D)
  if numel(D) == 1
    D = D*speye(length(y));
  end
  D = @(x) D\x;
end

if isnumeric(A)
  A = @(x) A*x;
end

z = D(y);
x = z;
for i=1:n
  z = z - D(A(z));
  x = x + z;
end
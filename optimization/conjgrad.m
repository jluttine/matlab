% x = conjgrad(A,b,x,I)
% A must be symmetric positive definite

function x = conjgrad(A,b,x,I,varargin)

options = struct( ...
    'maxiter', 100, ...
    'tol',     1e-6);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

if nargin < 4 || isempty(x)
  x = zeros(size(b));
end

if nargin < 5 || isempty(I)
  I = true(size(b));
end

if isnumeric(A)
  A = @(x) A*x;
end

r = b-A(x);
r(~I) = 0;

p = r;
rsold=r(I)'*r(I);

%Ap = zeros(size(x));

for i=1:options.maxiter
  Ap = A(p);
  Ap(~I) = 0;
  
  alpha = rsold/(p(:)'*Ap(:));
  x(:) = x(:) + alpha*p(:);
  r(:) = r(:) - alpha*Ap(:);
  rsnew = r(:)'*r(:);
  if sqrt(rsnew) < options.tol
    break;
  end
  p(:) = r(:) + (rsnew/rsold)*p(:);
  rsold = rsnew;
end

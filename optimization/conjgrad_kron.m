% x = conjgrad_kron(A1,A2,b,x,I)
% A1 and A2 must be symmetric positive definite

function x = conjgrad_kron(A1,A2,b,x,I,opts)

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

r=b-kronprod(A1,A2,x);
r(~I) = 0;

p=r;
rsold=r(I)'*r(I);

Ap = zeros(size(x));

for i=1:opts.maxiter
  Ap(:,:) = A2*p;
  Ap(:,:) = Ap*A1';
  
  alpha=rsold/(p(I)'*Ap(I));
  x(I)=x(I)+alpha*p(I);
  r(I)=r(I)-alpha*Ap(I);
  rsnew=r(I)'*r(I);
  if sqrt(rsnew) < opts.tol
    break;
  end
  p(I)=r(I)+(rsnew/rsold)*p(I);
  rsold=rsnew;
end

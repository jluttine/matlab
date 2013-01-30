% x = conjgradmv(A, b, I, varargin)

% Last modified 2011-02-25
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = conjgradmv(A, x, I, varargin)
% To reduce memory usage: b=x

% Large variables in memory: r, x, p, Ap

% $$$ function [x1, x1_mean] = gaussian_rand_pcg(x2, K, K1, z1, z2, M, S, varargin)

options = struct( ...
    'maxiter', 1000, ...
    'tol',     1e-6, ...
    'verbose', false, ...
    'x0',      [], ...
    'debug',   false);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% error(nargchk(5,7,nargin,'struct'));

% $$$ disp('Keyboard in conjgradmv');
% $$$ keyboard()

if nargin < 3 || isempty(I)
  I = true(size(x));
end
Imv = ~I;

if isnumeric(A)
  A = @(x) A*x;
end

% X(I) = Y uses a lot of memory, silly matlab..

x(Imv) = 0;
norm_b = sqrt(x(:)'*x(:));

% r = b - A*x
r = x;
if isempty(options.x0)
  x = zeros(size(x));
else
  x = options.x0;
end
r = r-A(x);
r(Imv) = 0;

p = r;
rsold=r(:)'*r(:);

if options.debug
  errors = nan(options.maxiter,1);
end

for i=1:options.maxiter
  Ap = A(p);
  Ap(Imv) = 0;
  
  alpha = rsold/(p(:)'*Ap(:));
  x = x + alpha*p;
  r = r - alpha*Ap;
  rsnew = r(:)'*r(:);
  res = sqrt(rsnew)/norm_b;
  if options.debug
    errors(i) = res;
  end
  if res < options.tol
    if options.verbose
      fprintf('Conjugate gradient converged in %d iterations to relative residual %e\n', ...
              i, res);
    end
    break;
  elseif i==options.maxiter && options.verbose
    fprintf(['Conjugate gradient did not converge in %d iterations having ' ...
             'relative residual %e\n'], options.maxiter, res);
  end
  p = r + (rsnew/rsold)*p;
  rsold = rsnew;
end
x(Imv) = 0;

if options.debug
  figure
  mean_error = mean(errors)
  semilogy(errors)
end

return

% $$$ 
% $$$ if isempty(S)
% $$$   x1 = z1; % L1(z1);
% $$$   z2 = x1 + z2 - x2 + mu2;
% $$$   %SKS = @(x) ( S*K(S'*x) );
% $$$   z2 = pcg(K, z2, options.tol, options.maxiter, M);
% $$$   x1 = mu1 + x1 - K1(z2);
% $$$   if nargout >= 2
% $$$     x1_mean = mu1 - K1( pcg(K, mu2-x2, options.tol, options.maxiter, M) );
% $$$   end
% $$$ elseif all(size(S)==size(x2)) && islogical(S)
% $$$   % S is logical array of missing values
% $$$   if ~isempty(M)
% $$$     M = @(x) remove(M(add(x,S)),S);
% $$$   end
% $$$   SKS = @(x) remove(K(add(x,S)),S);
% $$$   z2 = z1 + z2 - x2 + mu2;
% $$$   z2 = pcg(SKS, ...
% $$$           remove(z2, S), ...
% $$$           options.tol, ...
% $$$           options.maxiter, ...
% $$$           M);
% $$$   x1 = mu1 + z1 - K1(add(z2,S));
% $$$   if nargout >= 2
% $$$     x1_mean = mu1 - K1(add(pcg(SKS, ...
% $$$                                remove(mu2-x2,S), ...
% $$$                                options.tol, ...
% $$$                                options.maxiter, M), ...
% $$$                            S));
% $$$   end
% $$$ % $$$   disp('keyboard in gaussian_rand_pcg');
% $$$ % $$$   keyboard()
% $$$ else
% $$$   if ~isempty(M)
% $$$     M = @(x) S*M(S'*x);
% $$$   end
% $$$   SKS = @(x) ( S*K(S'*x) );
% $$$   z = pcg(SKS, S*(z1 + z2 - x2 + mu2), options.tol, options.maxiter, M);
% $$$   x1 = mu1 + z1 - K1(S'*z);
% $$$ 
% $$$   if nargout >= 2
% $$$     x1_mean = mu1 - K1( S'*pcg(SKS, S*(mu2-x2), options.tol, options.maxiter, M) );
% $$$   end
% $$$ end
% $$$ 
% $$$ function y = add(x,ind)
% $$$ % This also automatically reshapes to correct size
% $$$ y = zeros(size(ind));
% $$$ y(ind) = x;
% $$$ 
% $$$ function x = remove(x,ind)
% $$$ % This also automatically vectorizes
% $$$ x = x(ind);
% $$$ 
% $$$ function x = conjgrad(A,b,x,I,varargin)
% $$$ 
% $$$ options = struct( ...
% $$$     'maxiter', 100, ...
% $$$     'tol',     1e-6);
% $$$ 
% $$$ % Parse arguments
% $$$ [options, errmsg] = argparse( options, varargin{:} );
% $$$ error(errmsg);
% $$$ 
% $$$ if nargin < 4 || isempty(x)
% $$$   x = zeros(size(b));
% $$$ end
% $$$ 
% $$$ if nargin < 5 || isempty(I)
% $$$   I = true(size(b));
% $$$ end
% $$$ 
% $$$ if isnumeric(A)
% $$$   A = @(x) A*x;
% $$$ end
% $$$ 
% $$$ r = b-A(x);
% $$$ r(~I) = 0;
% $$$ 
% $$$ p = r;
% $$$ rsold=r(I)'*r(I);
% $$$ 
% $$$ %Ap = zeros(size(x));
% $$$ 
% $$$ for i=1:options.maxiter
% $$$   Ap = A(p);
% $$$   Ap(~I) = 0;
% $$$   
% $$$   alpha = rsold/(p(:)'*Ap(:));
% $$$   x(:) = x(:) + alpha*p(:);
% $$$   r(:) = r(:) - alpha*Ap(:);
% $$$   rsnew = r(:)'*r(:);
% $$$   if sqrt(rsnew) < options.tol
% $$$     break;
% $$$   end
% $$$   p(:) = r(:) + (rsnew/rsold)*p(:);
% $$$   rsold = rsnew;
% $$$ end

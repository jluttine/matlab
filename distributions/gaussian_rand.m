
function X = gaussian_rand(mu, Cov, N, varargin)

options = struct('type', 'cov');
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

if nargin < 3
  N = size(mu,2);
end

D = size(mu,1);

Z = randn(D,N);
switch options.type
 case 'cov'
  U = chol(Cov);
  Z = U' * Z;
 case 'invcov'
  U = chol(Cov);
  Z = linsolve_triu(U, Z);
 case 'cov-lchol'
  L = Cov;
  Z = L * Z;
 case 'invcov-lchol'
  Z = linsolve_triu(Cov, Z, true);
 otherwise
  error('Unknown type')
end
X = bsxfun(@plus, Z, mu);
%
% DATA
%

%randn('state', 10)
%rand('state', 10)
function test_gp_kron_eigs()

if nargin == 0
  N1 = 200;
  N2 = 300;
end
if nargin < 3
  N_samples = 4;
end
if nargin >= 4
  randn('state', seed);
  rand('state', seed);
else
  randn('state', 10);
  rand('state', 10);
end
  

% Inputs
x1 = 1:N1;
x2 = 1:N2;

% Covariance models

covfuncs = cell(2,3);
thetas = cell(2,3);
K = cell(2,3);

% Distances on the two domains
d1 = sqrt(sq_dist(x1(1), x1));
d2 = sqrt(sq_dist(x2(1), x2));

% "Spatial" covariance functions
covfunc = gp_cov_pp(d1,1);
covfunc = gp_cov_toeplitz(covfunc);
covfuncs(1,:) = {covfunc};
thetas(1,:) = {5; 10; 20};

% "Temporal" covariance functions
covfunc = gp_cov_pp(d2,1);
covfunc = gp_cov_toeplitz(covfunc);
covfuncs(2,:) = {covfunc};
thetas(2,:) = {5; 10; 20};

% Compute covariance matrices
for i=1:2
  for j=1:3
    K{i,j} = feval(covfuncs{i,j}, thetas{i,j});
  end
end


A = @(x) vec(kronprod(K{1,1},K{2,1},reshape(x,N2,N1)) + ...
             kronprod(K{1,2},K{2,2},reshape(x,N2,N1)) + ...
             kronprod(K{1,3},K{2,3},reshape(x,N2,N1))) + 0.3*x;

e = eigs(A,N1*N2, [], 1000, 'LM', struct('issym', true))
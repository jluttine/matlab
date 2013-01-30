
function test_lds(seed)

if nargin < 1
  seed = randn();
end
rand('state', seed);
randn('state', seed);

% Dynamical model (noisy resonator)
D = 2;
w = 0.2;
A = [cos(w) sin(w)/w
     -w*sin(w) cos(w)];
% $$$ A = [1 0
% $$$      -0.5 0.5];
q = randn(D);
Q = eye(D);%q*q';

% Observation model
M = 4;
C = randn(M,2);
r = randn(M);
s2 = 1e2;
R = s2*eye(M);%1e2*r*r';

% Initial state
mu_x0 = zeros(D,1);
Cov_x0 = eye(D);
x = mvnrnd(mu_x0,Cov_x0)';
%Covx0 = Covx0 * 1e10;

% Simulate data
N = 100;
Y = zeros(M,N);
X = zeros(D,N);
for n=1:N
  x = mvnrnd(A*x, Q)';
  Y(:,n) = mvnrnd(C*x, R)';
  X(:,n) = x;
end

% Missing vales
Y(:,30:60) = NaN;

%tsplot(Y)
%tsplot(X)
%X_true = X;

X_module = ...
    factor_module_lds(D, N, ...
                      'init', struct('mu_A', A, ...
                                     'Cov_A', 1e-5*eye(D)));

W_module = ...
    factor_module_ard(D, M, ...
                      'init', struct('mu_X', [], ...
                                     'Cov_X', 1e-1*eye(D)));

noise_module = ...
    noise_module_isotropic(M, N, ...
                           'prior', struct('a_tau', 1e6, ...
                                           'b_tau', s2*1e6));

Q = vbfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 50, ...
         'rotate', false, ...
         'update_w', true);

% $$$ C'
% $$$ Q.W_struct.posterior.mu_X
% $$$ Q.W_struct.posterior.Cov_X

mu_X = Q.X_struct.posterior.mu_X;
Cov_X = Q.X_struct.posterior.Cov_X;
var_X = reshape(Cov_X, [D*D, N]);
var_X = var_X(1:(D+1):end,:);
% $$$ tserrorplot(mu_X', 2*sqrt(var_X'))
% $$$ addtsplot(X, 'r');

mu_W = Q.W_struct.posterior.mu_X;
Cov_W = Q.W_struct.posterior.Cov_X;
[mu_Y, var_Y] = factor_product(mu_W, Cov_W, mu_X, Cov_X);
tserrorplot(mu_Y', 2*sqrt(var_Y'));
addtsplot(C*X, 'r');
addtsplot(Y, 'b+');

function [mu_Y, var_Y] = factor_product(mu_W, Cov_W, mu_X, Cov_X)
[D,M] = size(mu_W);
[D,N] = size(mu_X);

mu_Y = mu_W'*mu_X;
WW = bsxfun(@times, reshape(mu_W, [D 1 M]), reshape(mu_W, [1 D M])) + Cov_W;
XX = bsxfun(@times, reshape(mu_X, [D 1 N]), reshape(mu_X, [1 D N])) + Cov_X;
YY = reshape(WW, [D*D M])' * reshape(XX, [D*D N]);
var_Y = YY - mu_Y.^2;


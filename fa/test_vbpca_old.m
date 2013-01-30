function Q = test_vbpca(seed)

if nargin >= 1
  rand('state', seed);
  randn('state', seed);
end

D = 4;
M = 100;
N = 200;

X = randn(D,N);
W = diag(linspace(1,0.1,D))*randn(D,M);

Y = W'*X + 4;
Yn = Y + 0.1*randn(M,N);

% Estimate D+1 components

% PCA module (one constant component for bias modelling)
prior.mu = [zeros(D,1); 1];
prior.CovX = diag([ones(D,1); 1e-8]);
X_module = factor_module_iid(D+1,N,'prior',prior); % PCA module

% ARD module
W_module = factor_module_ard(D+1,M); 

% Isotropic noise
noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, 'init', 100);

Q = vbfa(D+1, Yn, W_module, X_module, noise_module, ...
          'maxiter', 50, ...
          'rotate', true, ...
          'rotation_checkgrad', false, ...
          'rotation_maxiter', 30);

noise_std = Q.Tau(1)^(-0.5)

%tsplot(Q.X)
%tsplot(Q.W)

figure
plot(Q.loglikelihood)

recon_error = rmse(Y, Q.W'*Q.X)

tsplot(Q.X);
tsplot(Q.W);

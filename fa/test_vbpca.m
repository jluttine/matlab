function Q = test_vbpca2(seed)

if nargin >= 1
  rand('state', seed);
  randn('state', seed);
end

D = 4;
M = 100;
N = 200;

X = randn(D,N);
W = diag(linspace(1,0.1,D))*randn(D,M);

Y = W'*X + 0;
Yn = Y + 0.1*randn(M,N);

% Missing values
%Yn(10:20,:) = NaN;
%Yn(:,10:20) = NaN;
Yn(rand(size(Yn))<0.2) = NaN;

% Estimate D+1 components
Dh = D + 1;

% PCA module (one constant component for bias modelling)
%prior.mu = [zeros(D+1,1)];
%prior.CovX = diag([ones(D+1,1)]);
prior.mu = [1; zeros(Dh,1)];
prior.CovX = diag([1e-3; ones(Dh,1)]);
X_noise = noise_module2_isotropic(N, 1, ...
                                  'init', struct('tau', 100), ...
                                  'prior', struct('a_tau', 1e-3, ...
                                                  'b_tau', 1e-3));
X_module = factor_module2_iid(Dh+1, N,...
                              'prior', prior, ...
                              'noise_module', []);

% ARD module
W_noise = noise_module2_isotropic(M, 1, ...
                                  'init', struct('tau', 100), ...
                                  'prior', struct('a_tau', 1e-3, ...
                                                  'b_tau', 1e-3));
W_module = factor_module2_ard(Dh+1, M, ...
                              'noise_module', W_noise); 

% Isotropic noise
noise_module = noise_module2_fixed(M,N,1);
%noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, 'init', 100);

Q = vbfa2(Dh+1, Yn, W_module, X_module, noise_module, ...
          'maxiter', 100, ...
          'rotate', true, ...
          'rotation_checkgrad', false, ...
          'rotation_maxiter', 30, ...
          'rotation_show', false);

Q = vbfa2(Dh+1, Yn, W_module, X_module, noise_module, ...
          'maxiter', 3, ...
          'rotate', true, ...
          'rotation_checkgrad', false, ...
          'rotation_maxiter', 30, ...
          'rotation_show', false);

%noise_std = Q.Tau(1)^(-0.5)

%tsplot(Q.X)
%tsplot(Q.W)

figure
plot(Q.loglikelihood)

recon_error = rmse(Y, Q.W'*Q.X)

tsplot(Q.X);
tsplot(Q.W);

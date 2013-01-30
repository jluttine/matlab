
%
% A simple demo about using VB for factor analysis models with transformation speed-ups.
%

% Number of features (rows)
M = 10;
% Number of samples (columns)
N = 200;
% Number of components (dimensionality of the latent space)
D = 3;


%
% Create artificial data
%
X = randn(D,N);
W = diag(linspace(1,0.1,D))*randn(D,M);
Y = W'*X;
% Add noise
Yn = Y + 2*randn(M,N);
% Add missing values
Yn(rand(size(Yn))<0.4) = NaN;

%
% Construct the model
%

% The number of components to estimate (in real applications, you don't know D)
Dh = D + 1;

% Prior model for X. The first component is a bias (mean 1 and small variance),
% the others have zero mean and unit variance.
prior.mu = [1; zeros(Dh,1)];
prior.CovX = diag([1e-3; ones(Dh,1)]);
X_module = factor_module_iid(Dh+1, N,...
                             'prior', prior, ...
                             'noise_module', []);

% Prior model for W.
% We include the noise parameter in the W module: it
% corresponds to a non-factorizing approximation Q(W,tau)
W_noise = noise_module_isotropic(M, 1, ...
                                 'init', struct('a_tau', 10, ...
                                                'b_tau', 1), ...
                                 'prior', struct('a_tau', 1e-3, ...
                                                 'b_tau', 1e-3));
% Use ARD for W in order to prune out unnecessary components
W_module = factor_module_ard(Dh+1, M, ...
                             'noise_module', W_noise); 
%W_module = factor_module_ard(Dh+1, M, ...
%                             'noise_module', []); 

% Prior model for isotropic noise.
% Fix this noise to 1 because the noise parameter is 
% actually included in the module for W.
noise_module = noise_module_fixed(M, N, 1);
% Alternatively, you could set noise_module of W to [] and then use
% the noise module below. That would correspond to a factorization 
% Q(X)*Q(W)*Q(tau).
%noise_module = noise_module_isotropic(M, N, ...
%                                      'init', struct('a_tau', 10, ...
%                                                     'b_tau', 1), ...
%                                      'prior', struct('a_tau', 1e-3, ...
%                                                      'b_tau', 1e-3));

%
% Run VB inference
%
% 'rotate' determines whether to use the speed-up transformations or not.
Q = vbfa(Dh+1, Yn, W_module, X_module, noise_module, ...
         'maxiter', 100, ...
         'rotate', true, ...
         'rotation_checkgrad', false, ...
         'rotation_maxiter', 30, ...
         'rotation_show', false);

%
% Show results
%

% Plot the VB lowerbound
figure
plot(Q.loglikelihood)
title('Log-likelihood lower bound')
xlabel('Iterations')

% Reconstruct
Yh = Q.W'*Q.X;
reconstruction_rmse = rmse(Y, Yh)

% Plot results
% Reconstruction
tsplot(Yh, 'k')
% True noiseless values
addtsplot(Y, 'b')
% Noisy observations
addtsplot(Yn, 'ro')

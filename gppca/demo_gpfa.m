
function demo_gpfa(seed)

if nargin >= 1
  randn('state', seed);
  rand('state', seed);
end

D = 4;

%
% Generate some data
%

%
% Generate X
%

covfunc_x = cell(D,1);

N_x = 300;
in_x = 1:N_x;

X = zeros(D,N_x);
D2_xx = sq_dist(in_x);
D_xx = sqrt(D2_xx);

pseudo_x = in_x(:,1:8:end);
D2_pp = sq_dist(pseudo_x);
D2_px = sq_dist(pseudo_x, in_x);
d2_x = diag(D2_xx);

is_pseudos_x = false(D,1);

% Trend component
covfunc = @(D2) gp_cov_se(D2);
covfunc_x{1} = gp_cov_pseudo(gp_cov_jitter(covfunc(D2_pp), 1e-3), ...
                             covfunc(D2_px), ...
                             covfunc(d2_x));
theta_x{1} = [100]; % true length scale
theta_init_x{1} = [50]; % init length scale
is_pseudos_x(1) = true;

% Almost periodic component (with fixed wave length)
covfunc = @(D2) gp_cov_product(gp_cov_periodic(sqrt(D2), ...
                                               'wavelength', 30), ...
                               gp_cov_se(D2));
covfunc_x{2} = gp_cov_pseudo(gp_cov_jitter(covfunc(D2_pp), 1e-3), ...
                             covfunc(D2_px), ...
                             covfunc(d2_x));
theta_x{2} = [1  % smoothness
              400]; % length scale of the decay
theta_init_x{2} = [0.5  % smoothness
                   800]; % length scale of the decay
is_pseudos_x(2) = true;

% Slow component
covfunc_x{3} = gp_cov_jitter(gp_cov_rq(D2_xx));
theta_x{3} = [5   % length scale
              1]; % alpha ("degrees of freedom")
theta_init_x{3} = [5     % length scale
                   1]; % alpha ("degrees of freedom")
is_pseudos_x(3) = false;

% Fast component
covfunc_x{4} = gp_cov_jitter(gp_cov_pp(D_xx, 1));
theta_x{4} = [4]; % length scale
theta_init_x{4} = [2.5]; % length scale
is_pseudos_x(4) = false;

% Generate X by random samples
for d=1:D
  if ~is_pseudos_x(d)
    K = covfunc_x{d}(theta_x{d});
    L = chol(K, 'lower');
    X(d,:) = L*randn(N_x,1);
  else
    [K_pp, K_px, k_x] = covfunc_x{d}(theta_x{d});
    L = chol(K_pp, 'lower');
    X(d,:) = K_px'*linsolve_lchol(L, L*randn(size(L,1),1));
  end    
end


%
% Generate W
%

N_w = 30;
in_w = rand(2,N_w);


D2_ww = sq_dist(in_w);
W = zeros(D,N_w);

% Meshgrid for showing the spatial functions
[grid_w1,grid_w2] = meshgrid(0:0.03:1, 0:0.03:1);
grid_w = [grid_w1(:)'; grid_w2(:)'];
D2_grid = sq_dist([grid_w, in_w]);
N_grid = size(grid_w, 2);
W_grid = zeros(D,N_grid);

% RQ covariance function for all spatial components
covfunc_w = cell(D,1);
covfunc_grid = cell(D,1);
covfunc_w(:) = {gp_cov_scale(gp_cov_jitter(gp_cov_rq(D2_ww)))};
covfunc_grid(:) = {gp_cov_scale(gp_cov_jitter(gp_cov_rq(D2_grid)))};
% True parameter values
theta_w{1} = [2     % signal magnitude
              0.4   % length scale
              3.0]; % alpha ("degrees of freedom")
theta_w{2} = [1     % signal magnitude
              0.2  % length scale
              2.0]; % alpha ("degrees of freedom")
theta_w{3} = [1     % signal magnitude
              0.1  % length scale
              1.0]; % alpha ("degrees of freedom")
theta_w{4} = [0.5   % signal magnitude
              0.05  % length scale
              1.0]; % alpha ("degrees of freedom")
% Parameter value initializations
theta_init_w = cell(D,1);
theta_init_w(:) = columns_to_cells(...
    [ones(1,D)            % signal magnitude
     linspace(0.3,0.01,D) % length scale
     ones(1,D)]);         % alpha ("degrees of freedom")
% $$$ theta_init_w(:) = {[1    % signal magnitude
% $$$                     0.02 % length scale
% $$$                     1]}; % alpha ("degrees of freedom")

% Generate W by random samples
for d=1:D
  K = covfunc_grid{d}(theta_w{d});
  L = chol(K, 'lower');
  w = L*randn(length(L),1);
  W_grid(d,:) = w(1:N_grid);
  W(d,:) = w((N_grid+1):end);
end

%
% Generate Y
%

% Noisy observations
Y_noiseless = W'*X;
s = 1;
Y = Y_noiseless + s*randn(N_w,N_x);

% Missing values
Y(rand(size(Y))<0.0) = nan;

%
% Plot  
%

% True spatial components
figure
for d=1:D
  subplot(ceil(sqrt(D)), ceil(sqrt(D)), d);
  % Plot contours of W
  cmax = max(abs(W_grid(d,:)));
  contourf(grid_w1, grid_w2, reshape(W_grid(d,:), size(grid_w1)));
  set(gca, 'clim', [-cmax, cmax])
  map_colormap();
  % Plot locations
  hold on
  plot(in_w(1,:), in_w(2,:), 'k+');
end

% True temporal components
tsplot(X);

%
% Inference with GPFA
%

% GP module with component-wise factorization for X
%theta_init_x = theta_x; % initialize with true parameter values
X_module = factor_module_gp_factorized(N_x, covfunc_x, theta_init_x, ...
                                       'update_hyperparameters', [5 10:10:100], ...
                                       'maxiter_hyperparameters', 10, ...
                                       'is_pseudo', is_pseudos_x, ...
                                       'init', zeros(D,N_x));

% GP module with component-wise factorization for W
%theta_init_w = theta_w; % initialize with true parameter values
W_module = factor_module_gp_factorized(N_w, covfunc_w, theta_init_w, ...
                                       'update_hyperparameters', [5 10:10:100], ...
                                       'maxiter_hyperparameters', 10);

% Isotropic noise
noise_module = noise_module_isotropic(N_w, N_x, ...
                                      'prior', struct('a_tau', 1e-3, ...
				                      'b_tau', 1e-3), ...
                                      'init', struct('a_tau', 10, ...
				                     'b_tau', 1));

% Run GPFA
Q = gpfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 50, ...
         'debug', false, ...
...%         'rotate', 1:10, ...
         'rotate', 1:15, ...
         'rotation_checkgrad', false, ...
         'rotation_show', false, ...
         'rotation_maxiter', 30);
% 'predict', true

fprintf('RMSE of the noiseless reconstruction: %f\n',  rmse(Y_noiseless, Q.W'*Q.X));
fprintf('STD of the noiseless predictions: %f\n', ...
        sqrt(mean(mean(Q.W.^2'*Q.CovX + Q.CovW'*Q.X.^2 + Q.CovW'*Q.CovX))));
fprintf('Estimated noise STD: %f\n', Q.Tau(1).^(-0.5))


% The development of the loglikelihood lower bound during the iteration
figure
plot(Q.loglikelihood)

% Estimated latent temporal signals
figure
for d=1:D
  subplot(D,1,d);
  errorplot(Q.X(d,:), 2*sqrt(Q.CovX(d,:)))
end



return



% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % $$$ tsplot(X)
% $$$ % $$$ tsplot(W)
% $$$ % $$$ tsplot(Y, 'k+')
% $$$ % $$$ return
% $$$ 
% $$$ % $$$ regularize = @(covfunc,N) gp_cov_sum(covfunc, ...
% $$$ % $$$                                      gp_cov_scale(gp_cov_delta(N), ...
% $$$ % $$$                                                   'scale', 1e-6));
% $$$ 
% $$$ %
% $$$ % GP module for X
% $$$ %
% $$$ 
% $$$ N_p = 2;
% $$$ pseudo_x = linspace(min(in_x),max(in_x), N_p);
% $$$ is_pseudo_x = false(D,1);
% $$$ 
% $$$ D2_xx = sq_dist(in_x);
% $$$ D2_pp = sq_dist(pseudo_x);
% $$$ D2_xp = sq_dist(in_x, pseudo_x);
% $$$ 
% $$$ % Covariance functions
% $$$ covfunc = @(D2) gp_cov_se(D2);
% $$$ covfunc_x = cell(D,1);
% $$$ covfunc_x{1} = gp_cov_jitter(gp_cov_se(D2_xx));
% $$$ % $$$ covfunc_x{1} = gp_cov_pseudo(gp_cov_jitter(gp_cov_se(D2_pp)), ...
% $$$ % $$$                              gp_cov_se(D2_xp), ...
% $$$ % $$$                              gp_cov_se(diag(D2_xx)));
% $$$ % $$$ is_pseudo_x(1) = true
% $$$ theta_x{1} = 30;
% $$$ 
% $$$ covfunc_x{2} = gp_cov_jitter(covfunc(D2_xx));
% $$$ theta_x{2} = 3;
% $$$ covfunc_x{3} = gp_cov_jitter(gp_cov_pp(sqrt(D2_xx),1));
% $$$ theta_x{3} = 1;
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ imagesc(covfunc_x{1}(theta_x{1}));
% $$$ % $$$ figure
% $$$ % $$$ imagesc(covfunc_x{2}(theta_x{2}));
% $$$ % $$$ return
% $$$ 
% $$$ % GP module with component-wise factorization
% $$$ X_module = factor_module_gp_factorized(N_x, ...
% $$$                                        covfunc_x, ...
% $$$                                        theta_x, ...
% $$$                                        'is_pseudo', is_pseudo_x, ...
% $$$                                        'update_hyperparameters', 2);
% $$$ 
% $$$ %
% $$$ % GP module for W
% $$$ %
% $$$ 
% $$$ D2_ww = sq_dist(in_w);
% $$$ 
% $$$ % Covariance functions
% $$$ covfunc_w = cell(D,1);
% $$$ covfunc_w{1} = gp_cov_scale(gp_cov_jitter(gp_cov_se(D2_ww)));
% $$$ theta_w{1} = [1; 4];
% $$$ covfunc_w{2} = gp_cov_scale(gp_cov_jitter(gp_cov_se(D2_ww)));
% $$$ theta_w{2} = [1; 3];
% $$$ covfunc_w{3} = gp_cov_scale(gp_cov_jitter(gp_cov_pp(sqrt(D2_ww),1)));
% $$$ theta_w{3} = [1; 2];
% $$$ 
% $$$ % GP module with component-wise factorization
% $$$ W_module = factor_module_gp_factorized(N_w, covfunc_w, theta_w, ...
% $$$                                        'update_hyperparameters', 2);
% $$$ 
% $$$ %
% $$$ % Isotropic noise module
% $$$ %
% $$$ 
% $$$ %noise_module = noise_module_fixed(1/s^2 * ones(N_w, N_x));
% $$$ noise_module = noise_module_isotropic(N_w, N_x, 1e-3, 1e-3, 'init', 100);
% $$$ 
% $$$ %
% $$$ % VB inference
% $$$ %
% $$$ 
% $$$ %
% $$$ % TODO:
% $$$ %
% $$$ % - pseudo inputs
% $$$ %
% $$$ % - put the noise to Q(W)
% $$$ %
% $$$ % - test more components
% $$$ %
% $$$ % - rotation, learn the hyperparameters jointly?
% $$$ %
% $$$ % - weighted noise
% $$$ %
% $$$ 
% $$$ Q = gpfa(D, Y, W_module, X_module, noise_module, ...
% $$$          'maxiter', 30, ...
% $$$          'update_noise', 2, ...
% $$$          'rotate', false);
% $$$ 
% $$$ noise_std = Q.Tau(1)^(-0.5)
% $$$ 
% $$$ tsplot(Q.X)
% $$$ tsplot(Q.W)
% $$$ %tsplot(Y, 'k+')
% $$$ 
% $$$ figure
% $$$ plot(Q.loglikelihood)
% $$$ 
% $$$ recon_error = rmse(Y_noiseless, Q.W'*Q.X)

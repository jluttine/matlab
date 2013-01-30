
function Q_gpfa = test_gpfa(seed)

if nargin >= 1
  randn('state', seed);
  rand('state', seed);
end

%
% Generate some data
%

D = 3;

in_x = 1:200;
N_x = length(in_x);
X = zeros(D,N_x);
X(1,:) = in_x/N_x;%sin(2*pi*in_w/10);
X(2,:) = cos(2*pi*in_x/7);
X(3,:) = randn(N_x,1);

in_w = 1:100;
N_w = length(in_w);
W = zeros(D,N_w);
W(1,:) = cos(2*pi*in_w/20);
W(2,:) = cos(2*pi*in_w/10);
W(3,:) = randn(N_w,1);

s = 0.5;
Y_noiseless = W'*X;
Y = Y_noiseless + s*randn(N_w,N_x);

Imv = (rand(size(Y)) < 0.4);
Imv(:,15) = true;
Y(Imv) = NaN;

%
% GP module for X
%

N_p = 10;
pseudo_x = linspace(min(in_x),max(in_x), N_p);
is_pseudo_x = false(D,1);

D2_xx = sq_dist(in_x);
D2_pp = sq_dist(pseudo_x);
D2_px = sq_dist(pseudo_x, in_x);

% Covariance functions
covfunc = @(D2) gp_cov_se(D2);
covfunc_x = cell(D,1);
% $$$ covfunc_x{1} = gp_cov_jitter(gp_cov_se(D2_xx));
covfunc_x{1} = gp_cov_pseudo(gp_cov_jitter(gp_cov_se(D2_pp)), ...
                             gp_cov_se(D2_px), ...
                             gp_cov_se(diag(D2_xx)));
is_pseudo_x(1) = true;
theta_x{1} = 30;

covfunc_x{2} = gp_cov_jitter(covfunc(D2_xx));
theta_x{2} = 3;

covfunc_x{3} = gp_cov_jitter(gp_cov_pp(sqrt(D2_xx),1));
theta_x{3} = 1.1;

% $$$ figure
% $$$ imagesc(covfunc_x{1}(theta_x{1}));
% $$$ figure
% $$$ imagesc(covfunc_x{2}(theta_x{2}));
% $$$ return


%
% GP module for W
%

D2_ww = sq_dist(in_w);

% Covariance functions
covfunc_w = cell(D,1);
covfunc_w{1} = gp_cov_scale(gp_cov_jitter(gp_cov_se(D2_ww)));
theta_w{1} = [1; 3];
covfunc_w{2} = gp_cov_scale(gp_cov_jitter(gp_cov_se(D2_ww)));
theta_w{2} = [1; 2];
covfunc_w{3} = gp_cov_scale(gp_cov_jitter(gp_cov_pp(sqrt(D2_ww),1)));
theta_w{3} = [1; 1.1];


%
% Isotropic noise module
%

%noise_module = noise_module_fixed(1/s^2 * ones(N_w, N_x));
noise_module = noise_module_isotropic(N_w, N_x, 1e-3, 1e-3, 'init', 100);

%
% VB inference
%

%
% TODO:
%
% - pseudo inputs
%
% - put the noise to Q(W)
%
% - test more components
%
% - rotation, optimize the hyperparameters jointly?
%
% - weighted noise
%

%
% GPFA
%

% GP modules with component-wise factorization

% NOTES:
%
% The algorithm can be quite sensitive to hyperparameter initialization and
% update schedule. Also, it seems that it's best to do rotations rarely and
% after the hyperparameters have been updated at least once. The rotations
% in GPFA are approximately optimized, so some problems can be caused by
% that. So it might be a good idea NOT to use rotations..

% Update only the smooth components at the beginning
maxiter = 30;
update_schedule = cell(maxiter,1);
update_schedule(:) = {1:D};
update_schedule(1:10) = {1};
%update_schedule(6:10) = {1:2};

X_module = factor_module_gp_factorized(N_x, covfunc_x, theta_x, ...
                                       'is_pseudo', is_pseudo_x, ...
                                       'update_hyperparameters', [5:5:1000], ...
                                       'update_schedule', {update_schedule});
W_module = factor_module_gp_factorized(N_w, covfunc_w, theta_w, ...
                                       'update_hyperparameters', [5:5:1000], ...
                                       'update_schedule', {update_schedule});

Q_gpfa = vbfa(D, Y, W_module, X_module, noise_module, ...
              'maxiter', maxiter, ...
              'update_noise', 1, ...
              'rotate', false, ...
              'rotation_checkgrad', false);


figure
plot(Q_gpfa.loglikelihood)

%
% PCA
%

X_module = factor_module_iid(D, N_x);
W_module = factor_module_ard(D, N_w);

Q_pca = vbfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 30, ...
         'update_noise', 1, ...
         'rotate', 1);


%
% Results
%

hax = tsplot(X);
title(hax(1), 'True X');
hax = tsplot(W);
title(hax(1), 'True W');

hax = tsplot(Q_gpfa.X);
title(hax(1), 'GPFA X');
hax = tsplot(Q_gpfa.W);
title(hax(1), 'GPFA W');

hax = tsplot(Q_pca.X);
title(hax(1), 'PCA X');
hax = tsplot(Q_pca.W);
title(hax(1), 'PCA W');

Yh_gpfa = Q_gpfa.W'*Q_gpfa.X;
Yh_pca = Q_pca.W'*Q_pca.X;
hax = tsplot([Y_noiseless(:,15)';
              Yh_gpfa(:,15)';
              Yh_pca(:,15)']);
title(hax(1), 'Reconstructions of observations');


% Comparison
%noise_std_gpfa = Q_gpfa.Tau(1)^(-0.5);
%noise_std_pca = Q_pca.Tau(1)^(-0.5);
recon_error_gpfa = rmse(Y_noiseless, Yh_gpfa)
recon_error_pca = rmse(Y_noiseless, Yh_pca)
test_error_gpfa = rmse(Y_noiseless(Imv), Yh_gpfa(Imv))
test_error_pca = rmse(Y_noiseless(Imv), Yh_pca(Imv))

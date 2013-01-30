function [Q, Yh] = hadsst3_experiment_gpfa_run()

% GP priors for loadings, PCA/isotropic prior for states.

% Number of components
D = 80;

%
% Process data
%

data = hadsst3_load_data();

% Form the data matrix
Y = data.observations;
[M,N] = size(Y);
Obs = ~isnan(Y);

%
% GP model for spatial W
%

ind = 0;

covfunc_w = cell(D,1);
theta_w = cell(D,1);
is_pseudos_w = false(D,1);

% Pseudo inputs (uniformly with respect to area size)
pseudo_w = points_on_sphere(18); % uniform points by number of latitudes
% Remove pseudo inputs that are on land (the nearest grid point is land)
ind_pseudo_w = mohsst5_points_to_grid_index(pseudo_w);
pseudo_w(:,mohsst5_is_land_index(ind_pseudo_w)) = [];
% $$$ % This code shows the pseudo inputs on the map
% $$$ figure
% $$$ map_projection('global-ellipse');
% $$$ map_plot(pseudo_w,'r+');
% $$$ map_coast()
% $$$ map_grid()
% $$$ return

% Transform inputs to 3-D Euclidean coordinates
in_w = data.coordinates;
in_w = geographic_to_euclidean(in_w);
pseudo_w = geographic_to_euclidean(pseudo_w);

% Squared distance matrices for the covariance functions
D2_ww = sq_dist(in_w);
D2_pp = sq_dist(pseudo_w);
D2_pw = sq_dist(pseudo_w, in_w);
d2_w = diag(D2_ww);

ind = 1:D;
fprintf('%d slow components for W (using %d pseudo inputs)\n', length(ind), ...
        size(pseudo_w,2));

% Covariance function (scaled squared exponential) with pseudo inputs
covfunc = @(D2) gp_cov_se(D2);
covfunc_w(ind) = {gp_cov_pseudo(...
    gp_cov_scale(gp_cov_jitter(covfunc(D2_pp), 1e-3)), ...
    gp_cov_scale(covfunc(D2_pw)), ...
    gp_cov_scale(covfunc(d2_w)))};

% Hyperparameters for the covariance functions
theta_w(ind) = columns_to_cells(...
    [linspace(1,0.1,length(ind));       % magnitudes
     linspace(5000,1000,length(ind))]); % lengthscales
is_pseudos_w(ind) = true;


%% Short scale components: piecewise polynomial in 3-D

ind = (ind(end)+1):D;
fprintf('%d fast components for W\n', length(ind));

% Use block-Toeplitz structure for the covariance function
[lat,lon0] = meshgrid(data.latitude,...
                      data.longitude(1));
in_w0 = geographic_to_euclidean([lon0(:)';lat(:)']);
d_ww = sqrt(sq_dist(in_w0, ...
                    geographic_to_euclidean(data.coordinates)));
covfunc = gp_cov_toeplitz_block(gp_cov_pp(d_ww,3));
% Add scaling and jitter
covfunc_w(ind) = {gp_cov_scale(gp_cov_jitter(covfunc))};
% Hyperparameters
theta_w(ind) = columns_to_cells(...
    [linspace(0.5,0.1,length(ind));     % magnitude
     linspace(3000,2000,length(ind))]); % lengthscale
is_pseudos_w(ind) = false;

%
% GPFA inference
%

% PCA module for X
%X_module = factor_module_iid(D, N);
covfunc_x = cell(D,1);
theta_x = cell(D,1);
covfunc_x(:) = {gp_cov_delta(N)};
theta_x(:) = {[]};
X_module = factor_module_gp_factorized(N, covfunc_x, theta_x, ...
                                       'update_hyperparameters', [], ...
                                       'init', zeros(D,N));

% Component-wise factorization for W
W_module = factor_module_gp_factorized(M, covfunc_w, theta_w, ...
                                       'update_hyperparameters', [5 10:10:100 100:25:2000], ...
                                       'maxiter_hyperparameters', 5, ...
                                       'is_pseudo', is_pseudos_w);

% Isotropic noise (precisions weighted proportionally to grid size)
weights = mohsst5_weights();
weights = repmat(weights, [1, N]);
noise_module = noise_module_product(...
    noise_module_isotropic(M, N, ...
                           'prior', struct('a_tau', 1e-3, ...
                                           'b_tau', 1e-3), ...
                           'init', struct('a_tau', 1, ...
                                          'b_tau', 0.01)), ...
    noise_module_fixed(M, N, weights));


% Filename for saving the results
folder = sprintf('/share/bayes/jluttine/results/hadsst3/gpfa');
mkdir(folder);
filename = sprintf('%s/results_hadsst3_gpfa_D=%d_%s', ...
                   folder, ...
                   D, ...
                   datestr(now,30));

% Run GPFA
Q = vbfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 200, ...
         'rotate', 1:50, ...
         'autosavefile', filename, ...
         'autosave', [1 5:10:2000]);

% Reconstruct
Yh = Q.W'*Q.X;

% Save the results
save(filename, '-struct', 'Q');

% Only short-scale GP components

function Q = mohsst5_experiment_gpfa(datanum)

%
% Load the data
%

[data,dataset,folder,maskfile] = metoffice_getdata(datanum, true);

% Number of components
D = 60;

%
% Model for temporal X
%

covfunc_x = cell(D,1);
theta_x = cell(D,1);
is_pseudos_x = false(D,1);


% Inputs (assume uniformly spaced time instances which is not exactly correct)
in_x = linspace(data.time(1), data.time(end), length(data.time));
% $$$ pseudo_x = in_x(1:10:end);

% $$$ % Squared distances for covariance functions
% $$$ D2_pp = sq_dist(pseudo_x);
% $$$ D2_px = sq_dist(pseudo_x, in_x);
% $$$ d2_x = zeros(size(in_x,2),1); %diag(D2_xx);
% $$$ D_pp = sqrt(D2_pp);
% $$$ D_px = sqrt(D2_px);
% $$$ d_x = sqrt(d2_x);

% Distance matrices for covariance functions
%d_xx = sqrt(D2_xx(1,:));
d_xx = sqrt(sq_dist(in_x(:,1),in_x));

% Fast components (4-18 months): piecewise polynomial in 1-D
% Take advantage of the Toeplitz structure of the covariance matrix
ind = 1:D;
fprintf('%d fast components for X\n', length(ind));
covfunc_x(ind) = {gp_cov_jitter(gp_cov_toeplitz(gp_cov_pp(d_xx,1)))};
theta_x(ind) = columns_to_cells(...
    [30*linspace(18,4,length(ind))]); % lengthscale or cut-off
is_pseudos_x(ind) = false;


%
% Model for spatial W
%

covfunc_w = cell(D,1);
theta_w = cell(D,1);
is_pseudos_w = false(D,1);

% Remove land area grid points
in_w = data.coordinates;
%in_w = mohsst5_remove_land(data.coordinates')';

%% Smooth components (using pseudo inputs)

% Pseudo inputs (uniformly with respect to area size)
pseudo_w = points_on_sphere(18); % uniform points by number of latitudes
% Remove pseudo inputs that are on land (the nearest grid point is land)
ind_pseudo_w = mohsst5_points_to_grid_index(pseudo_w);
mask = metoffice_get_mask(maskfile);
pseudo_w(:,~mask(ind_pseudo_w)) = [];

% $$$ % This code shows the pseudo inputs on the map
% $$$ figure
% $$$ map_projection('global-ellipse');
% $$$ map_plot(pseudo_w,'r+');
% $$$ map_coast()
% $$$ map_grid()
% $$$ return

% Transform inputs to 3-D Euclidean coordinates
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
                                 
% $$$ %% Short scale components: piecewise polynomial in 3-D
% $$$ 
% $$$ ind = (ind(end)+1):D;
% $$$ fprintf('%d fast components for W\n', length(ind));
% $$$ 
% $$$ % Use block-Toeplitz structure for the covariance function
% $$$ [lat,lon0] = meshgrid(data.latitude,data.longitude(1));
% $$$ in_w0 = geographic_to_euclidean([lon0(:)';lat(:)']);
% $$$ d_ww = sqrt(sq_dist(in_w0, ...
% $$$                     geographic_to_euclidean(data.coordinates)));
% $$$ covfunc = gp_cov_toeplitz_block(gp_cov_pp(d_ww,3));
% $$$ % Select only sea areas (discard land areas)
% $$$ sea = sum(~isnan(data.observations),2) > 0;
% $$$ covfunc = gp_cov_select(covfunc, sea);
% $$$ % Add scaling and jitter
% $$$ covfunc_w(ind) = {gp_cov_scale(gp_cov_jitter(covfunc))};
% $$$ % Hyperparameters
% $$$ theta_w(ind) = columns_to_cells(...
% $$$     [linspace(0.5,0.1,length(ind));     % magnitude
% $$$      linspace(3000,2000,length(ind))]); % lengthscale
% $$$ is_pseudos_w(ind) = false;

%
% Process data
%

% Form the data matrix
Y = data.data;
%Y = mohsst5_remove_land(data.data);
[M,N] = size(Y);
Obs = ~isnan(Y);


%
% GPFA inference
%

% Filename for saving the results
folder = [folder '/gpfa'];
mkdir(folder);
filename = sprintf('%s/results_rectest_%s_gpfa2_D=%d_anomalies=%d_%s', ...
                   folder, ...
                   dataset, ...
                   D, ...
                   anomalies, ...
                   datestr(now,'yyyymmdd'));

% Component-wise factorization for X
X_module = factor_module_gp_factorized(N, covfunc_x, theta_x, ...
                                       'update_hyperparameters', [5 10:10:100 100:25:2000], ...
                                       'maxiter_hyperparameters', 5, ...
                                       'is_pseudo', is_pseudos_x, ...
                                       'init', zeros(D,N));

% Component-wise factorization for W
W_module = factor_module_gp_factorized(M, covfunc_w, theta_w, ...
                                       'update_hyperparameters', [5 10:10:100 100:25:2000], ...
                                       'maxiter_hyperparameters', 5, ...
                                       'is_pseudo', is_pseudos_w);

% Isotropic noise (precisions weighted proportionally to grid size)
weights = repmat(data.gridsize, [1, N]);
% $$$ figure
% $$$ mohsst5_mapplot(metoffice_add_land(weights(:,1)));
% $$$ return
noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, ...
                                      'init', 10, ...
                                      'weights', weights);

% Run GPFA
Q = gpfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 2000, ...
         'rotate', 1:50, ... %[1:50 60:10:2000], ...
         'autosavefile', filename, ...
         'autosave', [10:100:2000]);

% Reconstruct
Yh = Q.W'*Q.X;

% Some performance measures
fprintf('Weighted training RMSE of the reconstruction: %f\n',  ...
        rmsew(Y(Obs)-Yh(Obs),weights(Obs)));

% Save the results
save(filename, '-struct', 'Q');

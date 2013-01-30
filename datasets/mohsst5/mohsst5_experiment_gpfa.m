function Q = mohsst5_experiment_gpfa()

% Number of components
D = 80;

data = mohsst5_loaddata();

%
% Model for temporal X
%

covfunc_x = cell(D,1);
theta_x = cell(D,1);
is_pseudos_x = false(D,1);


% Inputs (assume uniformly spaced time instances which is not exactly correct)
%in_x = data.time;
%pseudo_x = data.time(1:10:end);
in_x = linspace(data.time(1), data.time(end), length(data.time));
pseudo_x = in_x(1:10:end);

% Squared distances for covariance functions
%D2_xx = sq_dist(in_x);
D2_pp = sq_dist(pseudo_x);
D2_px = sq_dist(pseudo_x, in_x);
d2_x = zeros(size(in_x,2),1); %diag(D2_xx);

% Distance matrices for covariance functions
%d_xx = sqrt(D2_xx(1,:));
d_xx = sqrt(sq_dist(in_x(:,1),in_x));

% Slow components (1-20 years): rational quadratic using pseudo inputs
ind = 1:min(10,D);
fprintf('%d slow components for X\n', length(ind));
% $$$ covfunc = @(D2) gp_cov_se(D2); % SE
covfunc = @(D2) gp_cov_rq(D2); % RQ
covfunc_x(ind) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D2_pp), 1e-3), ...
                                covfunc(D2_px), ...
                                covfunc(d2_x))};
% $$$ theta_x(ind) = columns_to_cells(...
% $$$     [365*linspace(20,1,length(ind))]); % lengthscale for SE
theta_x(ind) = columns_to_cells(...
    [365*linspace(20,1,length(ind)) % lengthscale for RQ
     ones(1,length(ind))]);         % alpha for RQ
is_pseudos_x(ind) = true;

% Fast components (4-18 months): piecewise polynomial in 1-D
% Take advantage of the Toeplitz structure of the covariance matrix
ind = (ind(end)+1):D;
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
in_w = mohsst5_remove_land(data.coordinates')';

%% Smooth components (using pseudo inputs)

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
Y = mohsst5_remove_land(data.observations);
[M,N] = size(Y);
Obs = ~isnan(Y);

% Remove the test set from training data
Itest = load(sprintf('/share/bayes/data/jaakko/mohsst5/ind20test.mat'));
Itest = mohsst5_remove_land(Itest.Itest);
Itrain = ~Itest & Obs;
Itest = Itest & Obs;
Ytest = Y;
Ytest(~Itest) = nan;
Y(~Itrain) = nan;



%
% GPFA inference
%

% Filename for saving the results
folder = sprintf('/share/climate/jluttine/mohsst5/gpfa/%s', ...
                 datestr(now, 'yyyymmdd'));
mkdir(folder);
filename = sprintf('%s/results_mohsst5_gpfa_%s', ...
                   folder, ...
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
weights = mohsst5_weights();
weights = mohsst5_remove_land(weights);
weights = repmat(weights, [1, N]);
noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, ...
                                      'init', 10, ...
                                      'weights', weights);

% Run GPFA
Q = gpfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 2000, ...
         'rotate', 1:50, ...
         'autosavefile', filename, ...
         'autosave', [1 5:10:2000]);

% Reconstruct
Yh = Q.W'*Q.X;

% Some performance measures
sum_Itrain = sum(Itrain(:))
sum_Itest = sum(Itest(:))
fprintf('Weighted training RMSE of the reconstruction: %f\n',  ...
        rmsew(Y(Itrain)-Yh(Itrain),weights(Itrain)));
fprintf('Weighted test RMSE of the reconstruction: %f\n',  ...
        rmsew(Ytest(Itest)-Yh(Itest),weights(Itest)));
fprintf('STD of noisy predictions: %f\n', ...
        sqrt(mean(mean(Q.W.^2'*Q.CovX + Q.CovW'*Q.X.^2 + Q.CovW'*Q.CovX)) ...
             + mean(1./Q.Tau(:))));
%fprintf('Estimated noise STD: %f\n', Q.Tau(1).^(-0.5))

fprintf('Performance measure, weighted RMSE for the test set: %f\n', ...
        mohsst5_performance_rmsew(Yh, Itest));

% Save the results
save(filename, '-struct', 'Q');

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
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ 
% $$$ % Inputs
% $$$ in_x = data.time;
% $$$ pseudo_x = data.time(1:20:end);
% $$$ 
% $$$ % Squared distances for covariance functions
% $$$ D2_xx = sq_dist(in_x);
% $$$ D2_pp = sq_dist(pseudo_x);
% $$$ D2_xp = sq_dist(in_x, pseudo_x);
% $$$ d2_x = diag(D2_xx);
% $$$ 
% $$$ % Distance matrices for covariance functions
% $$$ D_xx = sqrt(D2_xx);
% $$$ D_pp = sqrt(D2_pp);
% $$$ D_xp = sqrt(D2_xp);
% $$$ d_x = sqrt(d2_x);
% $$$ 
% $$$ % Trend components: squared exponential
% $$$ covfunc = @(D2) gp_cov_se(D2);
% $$$ covfunc_x(1:5) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D2_pp), 1e-3), ...
% $$$                                 covfunc(D2_xp), ...
% $$$                                 covfunc(d2_x))};
% $$$ theta_x(1:5) = columns_to_cells(...
% $$$     [365*linspace(200,10,5)]); % lengthscale
% $$$ is_pseudos_x(1:5) = true;
% $$$ 
% $$$ 
% $$$ % Almost periodic components: periodic * squared exponential 
% $$$ % (fix period to one year)
% $$$ covfunc = @(D,D2) gp_cov_product(gp_cov_periodic(D, ...
% $$$                                                  'wavelength', 365), ...
% $$$                                  gp_cov_se(D2));
% $$$ covfunc_x(6:10) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D_pp,D2_pp), 1e-3), ...
% $$$                                  covfunc(D_xp,D2_xp), ...
% $$$                                  covfunc(d_x,d2_x))};
% $$$ theta_x(6:10) = columns_to_cells(...
% $$$     [ones(1,5);                 % smoothness
% $$$      365*linspace(300,100,5)]); % decay
% $$$ is_pseudos_x(6:10) = true;
% $$$ 
% $$$ % Smooth components: rational quadratic
% $$$ covfunc_x(11:15) = {gp_cov_jitter(gp_cov_rq(D2_xx))};
% $$$ theta_x(11:15) = columns_to_cells(...
% $$$     [365*linspace(10,1,5); % lengthscale
% $$$      ones(1,5)]);          % alpha
% $$$ is_pseudos_x(11:15) = false;
% $$$ 
% $$$ % Fast components: piecewise polynomial
% $$$ covfunc_x(16:D) = {gp_cov_jitter(gp_cov_pp(D_xx,1))};
% $$$ theta_x(16:D) = columns_to_cells(...
% $$$     [30*linspace(12,3,D-15)]); % lengthscale
% $$$ is_pseudos_x(16:D) = false;
% $$$ 
% $$$ 
% $$$ %
% $$$ % Model for spatial W
% $$$ %
% $$$ 
% $$$ covfunc_w = cell(D,1);
% $$$ theta_w = cell(D,1);
% $$$ is_pseudos_w = false(D,1);
% $$$ 
% $$$ % Locations
% $$$ in_w = geographic_to_euclidean(data.coordinates);
% $$$ 
% $$$ % Remove land areas
% $$$ land = sum(~isnan(data.observations),2) == 0;
% $$$ in_w(:,land) = [];
% $$$ 
% $$$ % Pseudo inputs
% $$$ % TODO: some elegant way to choose these!!!!!
% $$$ pseudo_w = in_w(:,1:10:end);
% $$$ 
% $$$ % Squared distance matrices
% $$$ D2_ww = sq_dist(in_w);
% $$$ D2_pp = sq_dist(pseudo_w);
% $$$ D2_wp = sq_dist(in_w, pseudo_w);
% $$$ d2_w = diag(D2_ww);
% $$$ 
% $$$ % Distance matrices
% $$$ D_ww = sqrt(D2_ww);
% $$$ 
% $$$ % Smooth components: rational quadratic (using pseudo inputs)
% $$$ covfunc = @(D2) gp_cov_rq(D2);
% $$$ covfunc_w(1:15) = {gp_cov_pseudo(...
% $$$     gp_cov_scale(gp_cov_jitter(covfunc(D2_pp), 1e-3)), ...
% $$$     gp_cov_scale(covfunc(D2_wp)), ...
% $$$     gp_cov_scale(covfunc(d2_w)))};
% $$$ theta_w(1:15) = columns_to_cells(...
% $$$     [linspace(1,0.1,15);      % magnitude
% $$$      linspace(10000,4000,15); % lengthscale
% $$$      ones(1,15)]);            % alpha
% $$$ is_pseudos_w(1:15) = true;
% $$$                                  
% $$$ % Short scale components: piecewise polynomial in 3-D
% $$$ covfunc = @(D) gp_cov_pp(D,3);
% $$$ covfunc_w(16:D) = {gp_cov_scale(gp_cov_jitter(covfunc(D_ww)))};
% $$$ theta_w(16:D) = columns_to_cells(...
% $$$     [linspace(0.5,0.1,D-15);     % magnitude
% $$$      linspace(6000,3000,D-15)]); % lengthscale
% $$$ is_pseudos_w(16:D) = false;
% $$$ 
% $$$ %
% $$$ % Process data
% $$$ %
% $$$ 
% $$$ Y = data.observations;
% $$$ Y(land,:) = [];
% $$$ [M,N] = size(Y);
% $$$ 
% $$$ %
% $$$ % GPFA inference
% $$$ %
% $$$ 
% $$$ % Component-wise factorization for X
% $$$ X_module = factor_module_gp_factorized(N, covfunc_x, theta_x, ...
% $$$                                        'update_hyperparameters', [1 5:5:100], ...
% $$$                                        'is_pseudo', is_pseudos_x);
% $$$ 
% $$$ % Component-wise factorization for W
% $$$ W_module = factor_module_gp_factorized(M, covfunc_w, theta_w, ...
% $$$                                        'update_hyperparameters', [1 5:5:100], ...
% $$$                                        'is_pseudo', is_pseudos_w);
% $$$ 
% $$$ % Isotropic noise
% $$$ noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, 'init', 10);
% $$$ 
% $$$ % Run GPFA
% $$$ Q = gpfa(D, Y, W_module, X_module, noise_module, ...
% $$$          'maxiter', 2, ...
% $$$          'update_noise', 1, ...
% $$$          'rotate', false, ...
% $$$          'rotation_checkgrad', false, ...
% $$$          'rotation_show', false, ...
% $$$          'debug', true);
% $$$ 
% $$$ fprintf('RMSE of the noiseless reconstruction: %f\n',  rmse(Y_noiseless, Q.W'*Q.X));
% $$$ fprintf('STD of the noiseless predictions: %f\n', ...
% $$$         sqrt(mean(mean(Q.W.^2'*Q.CovX + Q.CovW'*Q.X.^2 + Q.CovW'*Q.CovX))));
% $$$ fprintf('Estimated noise STD: %f\n', Q.Tau(1).^(-0.5))

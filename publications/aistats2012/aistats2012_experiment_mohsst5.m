
function [F,stdF,results] = aistats2012_experiment_mohsst5(model, testset, ...
                                                  seed)


% Seed for the random number generator
if nargin < 3
  seed = 10;
end
fprintf('Seed: %g\n', seed);
rand('state', seed);
randn('state', seed);

%
% FORM TRAIN AND TEST SETS
%

switch testset
 case 1 % randomly (uniformly) selected
  disp('Test set: Uniform');
  [Ytrain, Ytest, Itrain, Itest, data, sea] = aistats2012_mohsst5_sets(1,seed);
  testset_string = 'uniform';
 
 case 2 % pattern from earliest years used for latest years
  disp('Test set: Pattern');
  [Ytrain, Ytest, Itrain, Itest, data, sea] = aistats2012_mohsst5_sets(2,seed);
  testset_string = 'pattern';
  
 otherwise
  error('Unknown test set requested');
end

% Load the data
[M,N] = size(Ytrain);

switch model
 case 0
  model_string = 'gpkron';
 case 1
  model_string = 'gpkron-cov2';
 case 2
  model_string = 'vbpca';
 case 3
  model_string = 'gpfa';
 otherwise
  error('Unknown model requested');
end

filename = sprintf(['/home/jluttine/matlab/publications/aistats2012/' ...
                    'results_aistats2012_mohsst5_%s_%s_%s'], model_string, ...
                   testset_string, datestr(now,'yyyymmdd'));


%
% RUN ALGORITHM
%

switch model
 
 case {0,1}
 
  % 
  % Local GP with one or two covariance functions per domain
  %
  
  
  if model == 0
    disp('Model: Local GP with one covfuncs per domain');
    covfuncs_temporal = cell(1,1);
    covfuncs_spatial = cell(1,1);
  else
    disp('Model: Local GP with two covfuncs per domain');
    covfuncs_temporal = cell(1,2);
    covfuncs_spatial = cell(1,2);
  end
  % Two temporal covariance functions
  d = abs(1-(1:length(data.time)));
  covfuncs_temporal(:) = {gp_cov_toeplitz(gp_cov_pp(d,1))};

  % Two spatial covariance function
  [LON,LAT] = meshgrid(data.longitude,data.latitude);
  X = geographic_to_euclidean([LON(:)';LAT(:)']);
  % Use block-Toeplitz structure for the covariance function
  [lat,lon0] = meshgrid(data.latitude,data.longitude(1));
  X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
  D = sqrt(sq_dist(X0,X));
  covfunc = gp_cov_toeplitz_block(gp_cov_pp(D,3));
  % Select sea areas
  covfuncs_spatial(:) = {gp_cov_select(covfunc, sea)};

  if model == 0
    % Initial guess for covariance parameters
    theta_init = [0.5;  ... % signal 1 magnitude
                  10;    ... % signal 1 temporal length scale
                  3000; ... % signal 1 spatial length scale
                  0.5]';    % noise magnitude
  else
    % Initial guess for covariance parameters
    theta_init = [0.5;  ... % signal 1 magnitude
                  15;    ... % signal 1 temporal length scale
                  4000; ... % signal 1 spatial length scale
                  0.5;  ... % signal 2 magnitude
                  5;    ... % signal 2 temporal length scale
                  3000; ... % signal 2 spatial length scale
                  0.5]';    % noise magnitude
  end

  % Prior for the hyperparameters
  a = 1e-3;
  b = 1e-3;
  logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
  dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);
  
  % Noise scaling using area size
  w_spatial = 1./sqrt(cosd(LAT));
  %noise_scale = repmat(w(sea), [1,size(Y,2)]);

  % Inference
  N_samples = 1000;
  burnin = floor(N_samples/2);
  res = gp_kron(Ytrain, ...
                covfuncs_temporal, ...
                covfuncs_spatial, ...
                N_samples, ...
                theta_init, ...
                logprior_theta, ...
                dlogprior_theta, ...
                burnin, ...
                'noise_scale2', w_spatial(sea));

  % Results
  % F = sum(res.F,3);
  res
  if model == 0
    F = res.F;
    FF = res.FF;
  else
    F = res.sumF;
    FF = res.sumF2;
  end
  % Note, this doesn't really make any sense! You should compute like:
  % samplemean(sum(F,3).^2), not sum(samplemean(F.^2),3)
  stdF = sqrt(FF - F.^2);
  results = res;
  
  
 case 2
  
  %
  % VB PCA
  %
 
  disp('Model: VB PCA');
  
  D = 80;

  % PCA module for X (one constant component for modeling bias)
  prior.mu = [1; zeros(D-1,1)];
  prior.CovX = diag([1e-6; ones(D-1,1)]);
  X_module = factor_module_iid(D, N, 'prior', prior);

  % ARD module for W
  W_module = factor_module_ard(D, M);

  % Isotropic noise (precisions weighted proportionally to grid size)
  [LON,LAT] = meshgrid(data.longitude, data.latitude);
  weights = cosd(LAT(:));
  weights = weights(sea); %metoffice_remove_bins(weights,maskfile);
  weights = repmat(weights, [1, N]);
  noise_module = noise_module_product(...
      noise_module_isotropic(M, N, ...
                             'prior', struct('a_tau', 1e-3, ...
                                             'b_tau', 1e-3), ...
                             'init', struct('a_tau', 1e-3, ...
                                            'b_tau', 1e-3)), ...
      noise_module_fixed(M, N, weights));
% $$$   noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, ...
% $$$                                         'init', 10, ...
% $$$                                         'weights', weights);

  % Run VB PCA
  Q = vbfa(D, Ytrain, W_module, X_module, noise_module, ...
           'maxiter', 200, ...
           'rotate', true);
% $$$            'autosavefile', filename, ...
% $$$            'autosave', [1 20:20:2000]);

  % Results
  F = Q.W'*Q.X;
  varF = reshape(Q.CovW, [D*D, M])' * reshape(Q.CovX, [D*D,N]);
  for m=1:M
    for n=1:N
      varF(m,n) = varF(m,n) + ...
          Q.W(:,m)'*Q.CovX(:,:,n)*Q.W(:,m) + ...
          Q.X(:,n)'*Q.CovW(:,:,m)*Q.X(:,n);
    end
  end
  stdF = sqrt(varF);
  results = Q;

 case 3
  
  %
  % GPFA
  %
  

  % Number of components
  D = 80;

%  data = mohsst5_loaddata();

  %
  % Model for temporal X
  %

  covfunc_x = cell(D,1);
  theta_x = cell(D,1);
  is_pseudos_x = false(D,1);


  % Inputs (assume uniformly spaced time instances which is not exactly correct)
  in_x = linspace(data.time(1), data.time(end), length(data.time));
  pseudo_x = in_x(1:10:end);

  % Squared distances for covariance functions
  D2_pp = sq_dist(pseudo_x);
  D2_px = sq_dist(pseudo_x, in_x);
  d2_x = zeros(size(in_x,2),1); %diag(D2_xx);
  D_pp = sqrt(D2_pp);
  D_px = sqrt(D2_px);
  d_x = sqrt(d2_x);

  % Distance matrices for covariance functions
  %d_xx = sqrt(D2_xx(1,:));
  d_xx = sqrt(sq_dist(in_x(:,1),in_x));

  % Slow components (1-20 years): rational quadratic using pseudo inputs
  ind = 1:min(10,D);
  fprintf('%d slow components for X\n', length(ind));
  covfunc = @(D2) gp_cov_rq(D2); % RQ
  covfunc_x(ind) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D2_pp), 1e-6), ...
                                  covfunc(D2_px), ...
                                  covfunc(d2_x))};
  theta_x(ind) = columns_to_cells(...
      [365*linspace(5,1,length(ind)) % lengthscale for RQ
       ones(1,length(ind))]);         % alpha for RQ
  is_pseudos_x(ind) = true;
  
  % Quasi-periodic components
  ind = ind(end)+(1:5);
  fprintf('%d quasi-periodic components for X\n', length(ind));
  covfunc = @(D,D2) gp_cov_product(gp_cov_periodic(D, 'wavelength', 365), ...
                                   gp_cov_se(D2));
  covfunc_x(ind) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D_pp,D2_pp), 1e-6), ...
                                  covfunc(D_px,D2_px), ...
                                  covfunc(d_x,d2_x))};
  theta_x(ind) = columns_to_cells(...
      [ones(1,length(ind))           % smoothness of periodicity
       365*100*ones(1,length(ind))]); % lengthscale of SE decay
  is_pseudos_x(ind) = true;

  % Fast components (4-12 months): piecewise polynomial in 1-D
  % Take advantage of the Toeplitz structure of the covariance matrix
  ind = (ind(end)+1):D;
  fprintf('%d fast components for X\n', length(ind));
  covfunc_x(ind) = {gp_cov_jitter(gp_cov_toeplitz(gp_cov_pp(d_xx,1)))};
  theta_x(ind) = columns_to_cells(...
      [30*linspace(6,4,length(ind))]); % lengthscale or cut-off
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
  % uniform points by number of latitudes
  pseudo_w = points_on_sphere(18);
  % Remove pseudo inputs that are on land (the nearest grid point is land)
  ind_pseudo_w = mohsst5_points_to_grid_index(pseudo_w);
  pseudo_w(:,mohsst5_is_land_index(ind_pseudo_w)) = [];

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
  covfunc = @(D2) gp_cov_rq(D2);
  %covfunc = @(D2) gp_cov_se(D2);
  covfunc_w(ind) = {gp_cov_pseudo(...
      gp_cov_scale(gp_cov_jitter(covfunc(D2_pp), 1e-6)), ...
      gp_cov_scale(covfunc(D2_pw)), ...
      gp_cov_scale(covfunc(d2_w)))};

  % Hyperparameters for the covariance functions
  theta_w(ind) = columns_to_cells(...
      [linspace(1,0.1,length(ind));     % magnitudes
       linspace(3000,2000,length(ind)); % lengthscales
       linspace(3,3,length(ind))]);     % alpha for RQ
% $$$   theta_w(ind) = columns_to_cells(...
% $$$       [linspace(1,0.1,length(ind));       % magnitudes
% $$$        linspace(5000,1000,length(ind))]); % lengthscales
  is_pseudos_w(ind) = true;
  
  % GPFA inference

% $$$   % Filename for saving the results
% $$$   folder = sprintf('/share/climate/jluttine/mohsst5/gpfa/%s', ...
% $$$                    datestr(now, 'yyyymmdd'));
% $$$   mkdir(folder);
% $$$   filename = sprintf('%s/results_mohsst5_gpfa_%s', ...
% $$$                      folder, ...
% $$$                      datestr(now,'yyyymmdd'));

  % Component-wise factorization for X
  X_module = factor_module_gp_factorized(N, covfunc_x, theta_x, ...
                                         'update_hyperparameters', [3 5 10:10:100 100:25:2000], ...
                                         'maxiter_hyperparameters', 5, ...
                                         'is_pseudo', is_pseudos_x, ...
                                         'init', struct('X', zeros(D,N)));

  % Component-wise factorization for W
  W_module = factor_module_gp_factorized(M, covfunc_w, theta_w, ...
                                         'update_hyperparameters', [3 5 10:10:100 100:25:2000], ...
                                         'maxiter_hyperparameters', 5, ...
                                         'is_pseudo', is_pseudos_w);

  % Isotropic noise (precisions weighted proportionally to grid size)
  weights = mohsst5_weights();
  weights = mohsst5_remove_land(weights);
  weights = repmat(weights, [1, N]);
  noise_module = noise_module_product(...
      noise_module_isotropic(M, N, ...
                             'prior', struct('a_tau', 1e-3, ...
                                             'b_tau', 1e-3), ...
                             'init', struct('a_tau', 1e-3, ...
                                            'b_tau', 1e-3)), ...
      noise_module_fixed(M, N, weights));

  % Run GPFA
  Q = gpfa(D, Ytrain, W_module, X_module, noise_module, ...
           'maxiter', 200, ...
           'rotate', 1:50);
% $$$            'autosavefile', filename, ...
% $$$            'autosave', [1 5:10:2000]);

  % Reconstruct
  % Results
  F = Q.W'*Q.X;
  varF = Q.CovW' * Q.CovX;
  for m=1:M
    for n=1:N
      varF(m,n) = varF(m,n) + ...
          Q.W(:,m)'*diag(Q.CovX(:,n))*Q.W(:,m) + ...
          Q.X(:,n)'*diag(Q.CovW(:,m))*Q.X(:,n);
    end
  end
  stdF = sqrt(varF);
  results = Q;

 otherwise
  error('Unknown model requested');
end


%
% SHOW RESULTS
%

% Error measures
% $$$ [LON,LAT] = meshgrid(data.longitude, data.latitude);
% $$$ weights = cosd(LAT(:));
% $$$ weights = weights(sea);
% $$$ weights = repmat(weights, [1, N]);
% $$$ rmse_train = rmsew(Y(Itrain)-F(Itrain),weights(Itrain));
% $$$ rmse_test = rmsew(Y(Itest)-F(Itest),weights(Itest));
% $$$ rmse_zero = rmsew(Y(Itest),weights(Itest));
[rmse_train, rmse_test] = aistats2012_mohsst5_rmse(F, testset, seed);
fprintf('Training WRMSE=%.4f and testing WRMSE=%.4f\n', rmse_train, rmse_test);
% fprintf('Testing WRMSE=%.4f with zero predictions\n', rmse_zero);


save(filename, 'F', 'stdF', 'results', 'Ytest', 'Ytrain');
fprintf('Saved results to %s\n', filename);

if nargout < 1
  clear F
  clear stdF
  clear results
end

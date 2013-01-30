function F = metoffice_experiment_gpfa_gpkron(datanum, anomalies, validation)

%
% Load the data
%

date = datestr(now,'yyyymmdd');

if nargin < 3
  validation = true;
end

[data,dataset,folder,maskfile] = metoffice_getdata(datanum, anomalies, validation);

if anomalies
  disp('Model anomalies')
  comps = [0 1 0]; 
  comps_spatial = [0 1 0];
else
  disp('Don''t remove climatological averages')
  comps = [5 1 0]; 
  comps_spatial = [5 1 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPFA inference
%

debug = false;

% Number of components
D = sum(comps);
% DEBUGGING!!
if debug
  maxiter = 1;
  N_samples = 2;
else
  maxiter = 1000;
  N_samples = 2000;
end

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
d_xx = sqrt(sq_dist(in_x(:,1),in_x));

ind = 0;

% Periodic components (1 year period) with decay (rational quadratic)
% $$$ ind = 1:min(D,comps(1));
% $$$ fprintf('%d periodic components for X\n', length(ind));
% $$$ covfunc = @(D,D2) gp_cov_product(gp_cov_periodic(D, 'wavelength', 365), ...
% $$$                                  gp_cov_rq(D2));
% $$$ covfunc_x(ind) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D_pp,D2_pp), 1e-6), ...
% $$$                                 covfunc(D_px,D2_px), ...
% $$$                                 covfunc(d_x,d2_x))};
% $$$ theta_x(ind) = columns_to_cells(...
% $$$     [linspace(1,1,length(ind))      % smoothness of the period
% $$$      365*linspace(20,1,length(ind)) % lengthscale of the decay (RQ)
% $$$      ones(1,length(ind))]);         % alpha (RQ)
% $$$ is_pseudos_x(ind) = true;
% Periodic components (1 year period) WITHOUT decay
if comps(1) > 0
  ind = 1:min(D,comps(1));
  fprintf('%d periodic components for X\n', length(ind));
  covfunc = @(D) gp_cov_periodic(D, 'wavelength', 365.26);
  covfunc_x(ind) = {gp_cov_pseudo(gp_cov_jitter(covfunc(D_pp), 1e-6), ...
                                  covfunc(D_px), ...
                                  covfunc(d_x))};
  theta_x(ind) = columns_to_cells(...
      [linspace(1,1,length(ind))]);      % smoothness of the period
  is_pseudos_x(ind) = true;
end

% Slow components (1-20 years): rational quadratic using pseudo inputs
if comps(2) > 0
  ind = ind(end) + (1:min(D-ind(end),comps(2)));
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
end

% Fast components (4-18 months): piecewise polynomial in 1-D
% Take advantage of the Toeplitz structure of the covariance matrix
if comps(3) > 0
  ind = (ind(end)+1):D;
  fprintf('%d fast components for X\n', length(ind));
  covfunc_x(ind) = {gp_cov_jitter(gp_cov_toeplitz(gp_cov_pp(d_xx,1)))};
  theta_x(ind) = columns_to_cells(...
      [30*linspace(18,4,length(ind))]); % lengthscale or cut-off
  is_pseudos_x(ind) = false;
end

%
% Model for spatial W
%

covfunc_w = cell(D,1);
theta_w = cell(D,1);
is_pseudos_w = false(D,1);

% Remove land area grid points
in_w = data.coordinates;

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

ind = 0;

if comps_spatial(1) > 0
  ind = ind(end) + (1:comps_spatial(1));
  fprintf('%d iid components for W\n', length(ind));
  
  covfunc_w(ind) = {gp_cov_scale(gp_cov_delta(size(in_w,2)))};
  theta_w(ind) = columns_to_cells(...
      [linspace(1,0.1,length(ind))]);       % magnitudes
  
end

if comps_spatial(2) > 0
  ind = ind(end) + (1:comps_spatial(2));
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
       linspace(4000,1000,length(ind))]); % lengthscales
  is_pseudos_w(ind) = true;
end
%
% Process data
%

% Form the data matrix
Y = data.data;
[M,N] = size(Y);
Obs = ~isnan(Y);

% Filename for saving the results
folder = [folder '/gpfa_gpkron'];
mkdir(folder);
filename = sprintf('%s/results_rectest_%s_gpfa_D=%d_anomalies=%d_remval=%d_%s', ...
                   folder, ...
                   dataset, ...
                   D, ...
                   anomalies, ...
                   validation, ...
                   date);

% Component-wise factorization for X
X_module = factor_module_gp_factorized(N, covfunc_x, theta_x, ...
                                       'update_hyperparameters', [5 10:10:100 100:100:2000], ...
                                       'maxiter_hyperparameters', 5, ...
                                       'is_pseudo', is_pseudos_x, ...
                                       'init', zeros(D,N));

% Component-wise factorization for W
W_module = factor_module_gp_factorized(M, covfunc_w, theta_w, ...
                                       'update_hyperparameters', [5 10:10:100 100:100:2000], ...
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
         'maxiter', maxiter, ...
         'rotate', 1:100, ... %[1:50 60:10:2000], ...
         'autosavefile', filename, ...
         'autosave', [10:100:2000]);

Yh = Q.W'*Q.X;

% Some performance measures
fprintf('Weighted training RMSE of the reconstruction: %f\n',  ...
        rmsew(Y(Obs)-Yh(Obs),weights(Obs)));

% Save the results
save(filename, '-struct', 'Q');
fprintf('Saved GPFA results to %s\n', filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short-scale GP inference for the residuals
%

% Reconstruct
if debug
  Yh = zeros(size(Y));
else
  Yh = Q.W'*Q.X;
end
Yres = Y - Yh;
clear Q;


%
% Temporal covariance function (assume uniformly spaced time instances)
%

d = abs(1-(1:length(data.time)));
covfunc1 = gp_cov_toeplitz(gp_cov_pp(d,1));
theta_temporal = [7];
% $$$ covfunc1 = gp_cov_toeplitz(gp_cov_sum(gp_cov_pp(d,1), ...
% $$$                                       gp_cov_scale(gp_cov_pp(d,1))));
% $$$ theta_temporal = [7;   % length scale 1
% $$$                   1.0; % magnitude 2
% $$$                   3];  % length scale 2

%
% Spatial covariance function
%

[LON,LAT] = meshgrid(data.lon,data.lat);

X = geographic_to_euclidean([LON(:)';LAT(:)']);

% Use block-Toeplitz structure for the covariance function
[lat,lon0] = meshgrid(data.lat,data.lon(1));
X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
D = sqrt(sq_dist(X0,X));
covfunc2 = gp_cov_toeplitz_block(gp_cov_pp(D,3));
if datanum <= 4
  % 5x5 SST
  theta_spatial = [3000];  % length scale
elseif datanum == 5
  % 1x1 sea ice
  theta_spatial = [200];
end

% Select sea areas
sea = metoffice_get_mask(maskfile);
covfunc2 = gp_cov_select(covfunc2, sea);



%
% Inference
%

burnin = floor(N_samples/2);
folder_samples = sprintf('%s/samples_rectest_%s_gpkron_%s', ...
                         folder, ...
                         dataset, ...
                         date);
mkdir(folder_samples);
filename = sprintf('%s/results_rectest_%s_gpkron_D=%d_anomalies=%d_remval=%d_%s', ...
                   folder, ...
                   dataset, ...
                   sum(comps), ...
                   anomalies, ...
                   validation, ...
                   date);
filename_samples = sprintf('%s/samples_rectest_%s_gpkron_D=%d_anomalies=%d_remval=%d_%s', ...
                           folder_samples, ...
                           dataset, ...
                           sum(comps), ...
                           anomalies, ...
                           validation, ...
                           date);

a = 1e-3;
b = 1e-3;
logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);

% Initial guess for covariance parameters
theta_init = [0.5; ...               % total magnitude
              theta_temporal(:); ... % temporal parameters
              theta_spatial(:); ...  % spatial parameters
              0.5]';                 % noise magnitude

samplefunc = get_sample_function2(numel(theta_init), N_samples, burnin, ...
                                                filename, filename_samples);

% Weights for the noise levels using the respective grid size
w = 1./sqrt(cosd(LAT));
W = repmat(w(sea), [1,size(Y,2)]);

[get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron4(covfunc1, ...
                  solver_ldlchol(), ...
                  covfunc2, ...
                  solver_ldlchol(), ...
                  logprior_theta, ...
                  dlogprior_theta, ...
                  'samplefunc', samplefunc, ...
                  'rand_y', 'pcg', ...
                  'noise_scale', W, ...
                  'likelihood', 'whitened_prior');

% Transform to log-scale
func_theta = @(logtheta, varargin) func_theta_transformed(logtheta, ...
                                                  exp(logtheta), ...
                                                  func_theta, ...
                                                  varargin{:});
get_logpdf = @(f_theta) get_logpdf_transformed(f_theta, ...
                                               get_logpdf, ...
                                               sum(f_theta.theta_transformed));
get_dlogpdf = @(df_theta) get_dlogpdf_transformed(df_theta, ...
                                                  get_dlogpdf, ...
                                                  diag(exp(df_theta.theta_transformed)), ...
                                                  ones(size(df_theta.theta_transformed)));
theta_init = log(theta_init);

[rand_theta, f_theta] = mcmc_init_slicesampling(theta_init, ...
                                                get_logpdf, ...
                                                'fx', func_theta);

rand_model = get_rand_model(rand_theta, f_theta);




% Gibbs sampling
tic
gibbs(Yres,~Obs, N_samples, rand_model);
toc

% Check the results
res = samplefunc();

save(filename, '-struct', 'res');
disp(['Saved short-scale GP results to ', filename]);

% Mean reconstruction
F = Yh + res.F;

filename = sprintf('%s/results_rectest_%s_gpfa_gpkron_D=%d_anomalies=%d_remval=%d_%s', ...
                   folder, ...
                   dataset, ...
                   sum(comps), ...
                   anomalies, ...
                   validation, ...
                   date);
save(filename, 'F');
disp(['Saved total reconstruction to ', filename]);

fprintf('Weighted training RMSE of the reconstruction: %f\n',  ...
        rmsew(Y(Obs)-F(Obs),weights(Obs)));

if nargout < 1
  clear F;
end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = solver_ldlchol()

S.decompose = @decompose;
S.linsolve = @linsolve;
S.logdet = @logdet;
S.inv = @inv;
S.squareroot = @squareroot;
  
  function LD = decompose(K)
  [LD,p] = ldlchol(K);
  end
  
  function x = linsolve(LD, y)
  x = linsolve_ldlchol(LD,y);
  end
  
  function ldet = logdet(LD)
  ldet = logdet_ldlchol(LD);
  end
  
  function A = inv(LD)
  A = spinv_ldlchol(LD);
  end
  
  function L = squareroot(LD)
  L = ldlchol2lchol(LD);
  end

end


function [f_theta, df_theta] = func_theta_transformed(theta_transformed, ...
                                                  theta, func_theta, varargin)
if nargout <= 1
  f_theta = func_theta(theta, varargin{:});
  f_theta.theta_transformed = theta_transformed;
else
  [f_theta, df_theta] = func_theta(theta, varargin{:});
  f_theta.theta_transformed = theta_transformed;
  df_theta.theta_transformed = theta_transformed;
end
end

function logpdf_y = get_logpdf_transformed(fy, get_logpdf, logjacobian)
logpdf = get_logpdf(fy);
logpdf_y = @logpdf_transformed;
  function lpdf = logpdf_transformed(varargin)
  lpdf = logpdf(varargin{:}) + logjacobian;
  end
end

function dlogpdf_y = get_dlogpdf_transformed(dfy, get_dlogpdf, Jacobian, ...
                                                  dlogjacobian)
dlogpdf = get_dlogpdf(dfy);
dlogpdf_y = @dlogpdf_transformed;
  function dlpdf = dlogpdf_transformed(varargin)
  dlpdf = dlogpdf(varargin{:});
  dlpdf = Jacobian*dlpdf + dlogjacobian;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samplefunc = get_sample_function2(D_theta, N, burnin, filename, filename_samples)
results.Y = [];
results.F = 0;
results.FF = 0;
results.theta = zeros(D_theta, N);
samplefunc = @process_sample;
n = 1;
  function res = process_sample(Y, F, theta)
  if nargin >= 1
    % Store results
    results.Y = Y;
    if true && n > burnin
      results.F = (F + (n-burnin-1)*results.F) / (n-burnin);
    end
    if true && n > burnin
      results.FF = (F.*F + (n-burnin-1)*results.FF) / (n-burnin);
    end
    results.theta(:,n) = theta(:);
    % Save results
    fprintf('Saving results to %s..', filename)
    save(filename, '-struct', 'results');
    save(sprintf('%s_F%d',filename_samples,n), 'F');
    fprintf(' done.\n')
    n = n + 1;
  end
  if nargout >= 1
    % Return results
    res = results;
  end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gibbs(Y, Imv, N, rand_model, samplefunc)

Y(Imv) = 0;

for n=1:N
  
  t = cputime();
  Y = rand_model(Y,Imv);
  dt = cputime() - t;
  fprintf('Iteration step %d done. (%f seconds)\n', n, dt)
  
end

end



% Short-scale GP experiment for MetOffice data

%function metoffice_experiment_gpkron(datanum, validation)
function res = metoffice_experiment_gpkron(Y, times, coordinates, N_samples, ...
                                           filename, filename_samples)

%N_samples = 2000;

%
% Load the data
%

% $$$ if nargin < 2
% $$$   validation = true;
% $$$ end
% $$$ 
% $$$ [data,dataset,folder,maskfile] = metoffice_getdata(datanum, true, validation);

% Form the data matrix
% $$$ Y = data.data;
[M,N] = size(Y);
Obs = ~isnan(Y);

% Covariance type:
% 1) (k_1 * k_2) + noise
% 2) ((k_11 + k_12) * (k_21 + k_22)) + noise
type = 1;

%
% Temporal covariance function (assume uniformly spaced time instances)
%

d = abs(1-(1:length(times)));
switch type
 case 1
  covfunc1 = gp_cov_toeplitz(gp_cov_pp(d,1));
  theta_temporal = [7];   % length scale 1
 case 2
  covfunc1 = gp_cov_toeplitz(gp_cov_sum(gp_cov_pp(d,1), ...
                                        gp_cov_scale(gp_cov_pp(d,1))));
  theta_temporal = [7;   % length scale 1
                    1.0; % magnitude 2
                    3];  % length scale 2
end

%
% Spatial covariance function
%

%[LON,LAT] = meshgrid(data.lon,data.lat);
%X = geographic_to_euclidean([LON(:)';LAT(:)']);
X = geographic_to_euclidean(coordinates);

D = sqrt(sq_dist(X));
%%% Use block-Toeplitz structure for the covariance function
% $$$ [lat,lon0] = meshgrid(data.lat,data.lon(1));
% $$$ X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
% $$$ D = sqrt(sq_dist(X0,X));
switch type
 case 1
  covfunc2 = gp_cov_jitter(gp_cov_pp(D,3),1e-3);
  %covfunc2 = gp_cov_jitter(gp_cov_toeplitz_block(gp_cov_pp(D,3)),1e-3);
  %if datanum <= 4
    % 5x5 SST
    theta_spatial = [3000];  % length scale
% $$$   elseif datanum == 5
% $$$     % 1x1 sea ice
% $$$     theta_spatial = [200];
% $$$   end
 case 2
  covfunc2 = gp_cov_sum(gp_cov_pp(D,3), ...
                        gp_cov_scale(gp_cov_pp(D,3)));
% $$$   covfunc2 = gp_cov_toeplitz_block(gp_cov_sum(gp_cov_pp(D,3), ...
% $$$                                               gp_cov_scale(gp_cov_pp(D, ...
% $$$                                                     3))));
  covfunc2 = gp_cov_jitter(covfunc2, 1e-3);
  %if datanum <= 4
    % 5x5 SST
    theta_spatial = [3000;  % length scale 1
                     1.0;   % magnitude 2
                     2000]; % length scale 2
% $$$   elseif datanum == 5
% $$$     % 1x1 sea ice
% $$$     theta_spatial = [300;  % length scale 1
% $$$                      1.0;   % magnitude 2
% $$$                      100]; % length scale 2
% $$$   end
end

% Select sea areas
% $$$ sea = metoffice_get_mask(maskfile);
% $$$ covfunc2 = gp_cov_select(covfunc2, sea);



%
% Inference
%

burnin = floor(N_samples/2);
% $$$ folder = [folder '/gpkron'];
% $$$ folder_samples = sprintf('%s/samples_rectest_%s_gpkron_remval=%d_%s', ...
% $$$                          folder, ...
% $$$                          dataset, ...
% $$$                          validation, ...
% $$$                          datestr(now,'yyyymmdd'));
% $$$ mkdir(folder);
% $$$ mkdir(folder_samples);
% $$$ filename = sprintf('%s/results_rectest_%s_gpkron_remval=%d_%s', ...
% $$$                    folder, ...
% $$$                    dataset, ...
% $$$                    validation, ...
% $$$                    datestr(now,'yyyymmdd'));
% $$$ filename_samples = sprintf('%s/samples_rectest_%s_gpkron_remval=%d_%s', ...
% $$$                            folder_samples, ...
% $$$                            dataset, ...
% $$$                            validation, ...
% $$$                            datestr(now,'yyyymmdd'));

% Initial guess for covariance parameters
theta_init = [0.5; ...               % total magnitude
              theta_temporal(:); ... % temporal parameters
              theta_spatial(:); ...  % spatial parameters
              0.5]';                 % noise magnitude

a = 1e-3 * ones(size(theta_init));
b = 1e-3 * ones(size(theta_init));
% Prior for noise:
a(end) = 3;
b(end) = 6;
% Prior distributions
logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);

samplefunc = get_sample_function2(numel(theta_init), ...
                                  N_samples, ...
                                  burnin, ...
                                  filename, ...
                                  filename_samples);

% Weights for the noise levels using the respective grid size
%w = 1./sqrt(cosd(LAT));
%W = repmat(w(sea), [1,size(Y,2)]);
w = 1./sqrt(cosd(coordinates(2,:)));
%W = repmat(w(:), [1,size(Y,2)]);

res = gp_kron(Y, ...
              covfunc1, ...
              covfunc2, ...
              N_samples, ...
              theta_init, ...
              logprior_theta, ...
              dlogprior_theta, ...
              samplefunc, ...
              'noise_scale2', w)

% $$$ [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
% $$$     gp_init_kron4(covfunc1, ...
% $$$                   solver_ldlchol(), ...
% $$$                   covfunc2, ...
% $$$                   solver_ldlchol(), ...
% $$$                   logprior_theta, ...
% $$$                   dlogprior_theta, ...
% $$$                   'samplefunc', samplefunc, ...
% $$$                   'rand_y', 'pcg', ...
% $$$                   'noise_scale', W, ...
% $$$                   'likelihood', 'whitened_prior');

% $$$ % Transform to log-scale
% $$$ func_theta = @(logtheta, varargin) func_theta_transformed(logtheta, ...
% $$$                                                   exp(logtheta), ...
% $$$                                                   func_theta, ...
% $$$                                                   varargin{:});
% $$$ get_logpdf = @(f_theta) get_logpdf_transformed(f_theta, ...
% $$$                                                get_logpdf, ...
% $$$                                                sum(f_theta.theta_transformed));
% $$$ get_dlogpdf = @(df_theta) get_dlogpdf_transformed(df_theta, ...
% $$$                                                   get_dlogpdf, ...
% $$$                                                   diag(exp(df_theta.theta_transformed)), ...
% $$$                                                   ones(size(df_theta.theta_transformed)));
% $$$ theta_init = log(theta_init);
% $$$ 
% $$$ [rand_theta, f_theta] = mcmc_init_slicesampling(theta_init, ...
% $$$                                                 get_logpdf, ...
% $$$                                                 'fx', func_theta);
% $$$ 
% $$$ rand_model = get_rand_model(rand_theta, f_theta);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % Gibbs sampling
% $$$ tic
% $$$ gibbs(Y,~Obs, N_samples, rand_model);
% $$$ toc

% $$$ % Check the results
% $$$ res = samplefunc();

% $$$ time = 1600;
% $$$ figure
% $$$ 
% $$$ cl = max( max(abs(Y(:,time))), max(abs(res.F(:,time))) );
% $$$ clim = [-cl cl];
% $$$ 
% $$$ subplot(2,3,4)
% $$$ map_projection
% $$$ Z = nan(size(LON));
% $$$ Z(sea) = Y(:,time);
% $$$ map_pcolor(LON,LAT,Z);
% $$$ map_colormap()
% $$$ map_grid()
% $$$ map_coast()
% $$$ set(gca, 'clim', clim);
% $$$ 
% $$$ subplot(2,3,5)
% $$$ map_projection
% $$$ Z = nan(size(LON));
% $$$ Z(sea) = res.F(:,time);
% $$$ map_pcolor(LON,LAT,Z);
% $$$ map_colormap()
% $$$ map_grid()
% $$$ map_coast()
% $$$ set(gca, 'clim', clim);
% $$$ 
% $$$ subplot(2,3,6)
% $$$ map_projection
% $$$ std_F = sqrt(res.FF(:,time) - res.F(:,time).*res.F(:,time));
% $$$ Z = nan(size(LON));
% $$$ Z(sea) = std_F;
% $$$ map_pcolor(LON,LAT,Z);
% $$$ map_colormap()
% $$$ map_grid()
% $$$ map_coast()
% $$$ 
% $$$ subplot(2,3,1)
% $$$ semilogy(res.theta');
% $$$ subplot(2,3,2)
% $$$ plot(acorr(res.theta(:,(burnin+1):end)'));
% $$$ 
% $$$ location = 1000;
% $$$ subplot(2,3,3)
% $$$ plot([Y(location,:)', res.F(location,:)']);
% $$$ 
% $$$ % Theta scatter plot
% $$$ figure()
% $$$ plot_scatterhist(res.theta(:,(burnin+1):end)');
% $$$ 
% $$$ 
% $$$ % RMSE
% $$$ rmse_train = rmse(Y(Obs),res.F(Obs))
% $$$ rmse_test = rmse(Y(Obs),res.F(Obs))
% $$$ 
% $$$ theta = res.theta;
% $$$ save(filename, '-struct', 'res');
% $$$ disp(['Saved results to ', filename]);
% $$$ 

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ function S = solver_ldlchol()
% $$$ 
% $$$ S.decompose = @decompose;
% $$$ S.linsolve = @linsolve;
% $$$ S.logdet = @logdet;
% $$$ S.inv = @inv;
% $$$ S.squareroot = @squareroot;
% $$$   
% $$$   function S = decompose(K)
% $$$   [S.LD,p,S.q] = ldlchol(K);
% $$$   if p>0
% $$$     error('Matrix must be positive definite');
% $$$   end
% $$$   end
% $$$   
% $$$   function x = linsolve(S, y)
% $$$   x(S.q,:) = linsolve_ldlchol(S.LD,y(S.q,:));
% $$$   end
% $$$   
% $$$   function ldet = logdet(S)
% $$$   ldet = logdet_ldlchol(S.LD);
% $$$   end
% $$$   
% $$$   function A = inv(S)
% $$$   A(S.q,S.q) = spinv_ldlchol(S.LD);
% $$$   end
% $$$   
% $$$   function L = squareroot(S)
% $$$   L(S.q,S.q) = ldlchol2lchol(S.LD);
% $$$   end
% $$$ 
% $$$ end
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
    if ~isempty(filename)
      fprintf('Saving results to %s..', filename)
      save(filename, '-struct', 'results');
    end
    if ~isempty(filename_samples)
      save(sprintf('%s_F%d',filename_samples,n), 'F');
    end
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



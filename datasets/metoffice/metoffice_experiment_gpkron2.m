% Short-scale GP experiment for MetOffice data

function metoffice_experiment_gpkron2(datanum)

N_samples = 1000;

%
% Load the data
%

[data,dataset,folder,maskfile] = metoffice_getdata(datanum, true);

% Form the data matrix
Y = data.data;
[M,N] = size(Y);
Obs = ~isnan(Y);

% Covariance type:
% 1) ((k_11 + k_12) * (k_21 + k_22)) + noise
% 2) ((k_1 + noise) * (k_2 + noise)) + noise
type = 1;

%
% Temporal covariance function (assume uniformly spaced time instances)
%

d = abs(1-(1:length(data.time)));
covfunc1 = gp_cov_toeplitz(gp_cov_pp(d,1));
theta_temporal = [7]; % length scale

%
% Spatial covariance function
%

[LON,LAT] = meshgrid(data.lon,data.lat);

X = geographic_to_euclidean([LON(:)';LAT(:)']);

% Use block-Toeplitz structure for the covariance function
[lat,lon0] = meshgrid(data.lat,data.lon(1));
X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
D = sqrt(sq_dist(X0,X));
covfunc2 = gp_cov_jitter(gp_cov_toeplitz_block(gp_cov_pp(D,3)),1e-3);
% $$$ covfunc2 = gp_cov_toeplitz_block(gp_cov_sum(gp_cov_pp(D,3), ...
% $$$                                             gp_cov_scale(gp_cov_pp(D,3))));

% Select sea areas
sea = metoffice_get_mask(maskfile);
covfunc2 = gp_cov_select(covfunc2, sea);

if datanum <= 4
  % 5x5 SST
  theta_spatial = [3000];  % length scale
elseif datanum == 5
  % 1x1 sea ice
  theta_spatial = [50];
else
  error('Unknown datanum');
end


% $$$ K2 = covfunc2(theta_spatial);
% $$$ condest(K2)
% $$$ nnz(K2) / numel(K2)
% $$$ length(K2) / numel(K2)
% $$$ return

% $$$ % Re-order
% $$$ N_lon = length(data.lon);
% $$$ N_lat = length(data.lat);
% $$$ mask = metoffice_get_mask(maskfile);
% $$$ q = reshape(reshape(1:(N_lon*N_lat), [N_lon N_lat])', [N_lon*N_lat, 1]);
% $$$ [~,q] = sort(q(mask));
% $$$ K2 = covfunc2(theta_spatial);
% $$$ K2 = K2 + speye(size(K2));
% $$$ figure
% $$$ spy(K2)
% $$$ figure
% $$$ spy(K2(q,q));
% $$$ figure
% $$$ [L,p] = chol(K2);
% $$$ spy(L);
% $$$ figure
% $$$ [L,p,r] = chol(K2);
% $$$ spy(L);
% $$$ figure
% $$$ [L,p] = chol(K2(q,q));
% $$$ spy(L)
% $$$ %figure
% $$$ %imagesc(q'*K2*q)
% $$$ %figure
% $$$ %spy(q'*K2*q)
% $$$ return


%
% Inference
%

burnin = floor(N_samples/2);
folder = [folder '/gpkron2'];
folder_samples = sprintf('%s/samples_rectest_%s_gpkron2_%s', ...
                         folder, ...
                         dataset, ...
                         datestr(now,'yyyymmdd'));
mkdir(folder);
mkdir(folder_samples);
filename = sprintf('%s/results_rectest_%s_gpkron2_%s', ...
                   folder, ...
                   dataset, ...
                   datestr(now,'yyyymmdd'));
filename_samples = sprintf('%s/samples_rectest_%s_gpkron2_%s', ...
                           folder_samples, ...
                           dataset, ...
                           datestr(now,'yyyymmdd'));

a = 1e-3;
b = 1e-3;
logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);

% Initial guess for covariance parameters
theta_init = [0.5; ...               % total magnitude
              theta_temporal(:); ... % temporal parameters
              0.5; ...               % temporal noise magnitude
              theta_spatial(:); ...  % spatial parameters
              0.5]';                 % spatial noise magnitude

samplefunc = get_sample_function2(numel(theta_init), N_samples, burnin, ...
                                                filename, filename_samples);

% Weights for the noise levels using the respective grid size
w = 1./sqrt(cosd(LAT(:)));
w = w(sea);
%W = repmat(w(sea), [1,size(Y,2)]);

[get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron1(covfunc1, ...
                  solver_ldlchol(), ...
                  covfunc2, ...
                  solver_ldlchol(), ...
                  logprior_theta, ...
                  dlogprior_theta, ...
                  samplefunc, ...
                  'noise_scale2', w(:));

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
gibbs(Y,~Obs, N_samples, rand_model);
toc

% Check the results
res = samplefunc();

time = 1600;
figure

cl = max( max(abs(Y(:,time))), max(abs(res.F(:,time))) );
clim = [-cl cl];

subplot(2,3,4)
map_projection
Z = nan(size(LON));
Z(sea) = Y(:,time);
map_pcolor(LON,LAT,Z);
map_colormap()
map_grid()
map_coast()
set(gca, 'clim', clim);

subplot(2,3,5)
map_projection
Z = nan(size(LON));
Z(sea) = res.F(:,time);
map_pcolor(LON,LAT,Z);
map_colormap()
map_grid()
map_coast()
set(gca, 'clim', clim);

subplot(2,3,6)
map_projection
std_F = sqrt(res.FF(:,time) - res.F(:,time).*res.F(:,time));
Z = nan(size(LON));
Z(sea) = std_F;
map_pcolor(LON,LAT,Z);
map_colormap()
map_grid()
map_coast()

subplot(2,3,1)
semilogy(res.theta');
subplot(2,3,2)
plot(acorr(res.theta(:,(burnin+1):end)'));

location = 1000;
subplot(2,3,3)
plot([Y(location,:)', res.F(location,:)']);

% Theta scatter plot
figure()
plot_scatterhist(res.theta(:,(burnin+1):end)');


% RMSE
rmse_train = rmse(Y(Obs),res.F(Obs))
rmse_test = rmse(Y(Obs),res.F(Obs))

theta = res.theta;
save(filename, '-struct', 'res');
disp(['Saved results to ', filename]);


end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = solver_ldlchol()

S.decompose = @decompose;
S.linsolve = @linsolve;
S.logdet = @logdet;
S.inv = @inv;
S.squareroot = @squareroot;
  
  function S = decompose(K)
  [S.LD,p,S.q] = ldlchol(K);
  if p>0
    error('Matrix must be positive definite');
  end
  end
  
  function x = linsolve(S, y)
  x(S.q,:) = linsolve_ldlchol(S.LD,y(S.q,:));
  end
  
  function ldet = logdet(S)
  ldet = logdet_ldlchol(S.LD);
  end
  
  function A = inv(S)
  %A(S.q,S.q) = spinv_ldlchol(S.LD);
  N = length(S.q);
  R = sparse(1:N,S.q,ones(N,1));
  L = R'*spinv_ldlchol(S.LD)*R;
  end
  
  function L = squareroot(S)
  % WARNING/TODO: This permutation might take A LOT of time..
  % L(S.q,S.q) = ldlchol2lchol(S.LD);
  % Try instead:
  N = length(S.q);
  R = sparse(1:N,S.q,ones(N,1));
  L = R'*ldlchol2lchol(S.LD)*R;
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



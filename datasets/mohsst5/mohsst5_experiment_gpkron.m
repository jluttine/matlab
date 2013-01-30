
function mohsst5_experiment_gpkron(data)

if nargin < 1
  data = mohsst5_loaddata();
end

%
% Temporal covariance function
%

d = abs(1-(1:length(data.time)));
covfunc1 = gp_cov_toeplitz(gp_cov_sum(gp_cov_pp(d,1), ...
                                      gp_cov_scale(gp_cov_pp(d,1))));
theta_temporal = [7;   % length scale 1
                  1.0; % magnitude 2
                  3];  % length scale 2
%theta1 = 5;

%
% Spatial covariance function
%

[LON,LAT] = meshgrid(data.longitude,data.latitude);

X = geographic_to_euclidean([LON(:)';LAT(:)']);

% Use block-Toeplitz structure for the covariance function
[lat,lon0] = meshgrid(data.latitude,data.longitude(1));
X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
D = sqrt(sq_dist(X0,X));
covfunc2 = gp_cov_toeplitz_block(gp_cov_sum(gp_cov_pp(D,3), ...
                                            gp_cov_scale(gp_cov_pp(D,3))));

% Select sea areas
sea = sum(~isnan(data.observations),2) > 0;
covfunc2 = gp_cov_select(covfunc2, sea);

theta_spatial = [3000;  % length scale 1
                 1.0;   % magnitude 2
                 2000]; % length scale 2

%theta2 = 3000;

% DEBUG STUFF
K1 = covfunc1(theta_temporal);
K2 = covfunc2(theta_spatial);
whos
return

%
% Data
% 

Y = data.observations(sea,:);

% Divide data into test and train sets
testsize = 20; % testset size in percents
Itest = load(sprintf('/share/bayes/data/jaakko/mohsst5/ind%dtest.mat', ...
                     testsize));
Itest = Itest.Itest(sea,:) & ~isnan(Y);
Itrain = load(sprintf('/share/bayes/data/jaakko/mohsst5/ind%dtrain.mat', ...
                      testsize));
Itrain = Itrain.Itrain(sea,:) & ~isnan(Y);
Ytest = Y;
Ytest(Itrain) = nan;
Ytrain = Y;
Ytrain(Itest) = nan;


%
% Inference
%

N_samples = 2000;
burnin = floor(N_samples/2);
folder = sprintf('/share/climate/jluttine/mohsst5/gpkron/%s', ...
                 datestr(now, 'yyyymmdd'));
% $$$ folder = sprintf('/home/jluttine/matlab/icml2011/results/%s', datestr(now, ...
% $$$                                                   'yyyymmdd'));
mkdir(folder);
mkdir([folder '/samples']);
filename = sprintf('%s/results_mohsst5_gpkron_%s', folder, ...
                   datestr(now,'yyyymmdd'));
filename_samples = sprintf('%s/samples/samples_mohsst5_gpkron_%s', ...
                           folder, ...
                           datestr(now,'yyyymmdd'));
% $$$ burnin = min(100,floor(N_samples/2));

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
gibbs(Ytrain,~Itrain, N_samples, rand_model);
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
rmse_train = rmse(Y(Itrain),res.F(Itrain))
rmse_test = rmse(Y(Itest),res.F(Itest))
performance_train = mohsst5_performance_rmsew(res.F, Itrain)
performance_test = mohsst5_performance_rmsew(res.F, Itest)

% $$$ F = nan(size(data.observations));
% $$$ FF = nan(size(data.observations));
% $$$ Y = nan(size(data.observations));
% $$$ F(sea,:) = res.F;
% $$$ FF(sea,:) = res.FF;
% $$$ Y(sea,:) = res.Y;
theta = res.theta;
save(filename, '-struct', 'res');
disp(['Saved results to ', filename]);

% $$$ K = covfunc(theta_init);
% $$$ 
% $$$ figure
% $$$ k = zeros(size(LON));
% $$$ k(sea) = K(500,:);
% $$$ map_projection
% $$$ map_pcolor(LON,LAT,k);
% $$$ map_colormap()
% $$$ map_grid()
% $$$ map_coast()


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



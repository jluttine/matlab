
function [F,stdF,results] = nips2011_experiment_mohsst5(model, testset, seed)

% Seed for the random number generator
if nargin < 3
  seed = 10;
end
fprintf('Seed: %g\n', seed);
rand('state', seed);
randn('state', seed);

% Load the data
data = mohsst5_loaddata();
sea = sum(~isnan(data.observations),2) > 0;
Y = data.observations(sea,:);
[M,N] = size(Y);

switch testset
 case 1 % randomly (uniformly) selected
  disp('Test set: Uniform');
 
  I = rand(size(Y)) < 0.2;
  Itest = I & ~isnan(Y);
  Itrain = ~I & ~isnan(Y);
  Ytest = Y;
  Ytrain = Y;
  Ytest(~Itest) = nan;
  Ytrain(~Itrain) = nan;
  
  test_percent = sum(Itest(:)) / sum(~isnan(Y(:)));
  fprintf('Size of the test set: %.1f percent\n', 100*test_percent);

% $$$   figure
% $$$   spy(Itest)
% $$$   figure
% $$$   spy(Itrain)
% $$$   figure
% $$$   spy(~isnan(Y));
% $$$   figure
% $$$   spy(Itrain&Itest)
% $$$   figure
% $$$   spy((Itrain|Itest)==~isnan(Y));
% $$$   return
  
  testset_string = 'uniform';
 
 case 2 % pattern from earliest years used for latest years
  disp('Test set: Pattern');
  
  Imv = isnan(Y);
  Imv(:,(20*12+1):end) = false;
  Itest = Imv(:,end:-1:1) & ~isnan(Y);
  Itrain = ~Imv(:,end:-1:1) & ~isnan(Y);
  
  Ytest = Y;
  Ytrain = Y;
  Ytest(~Itest) = nan;
  Ytrain(~Itrain) = nan;
  
  test_percent = sum(Itest(:)) / sum(~isnan(Y(:)));
  fprintf('Size of the test set: %.1f percent\n', 100*test_percent);

% $$$   figure
% $$$   spy(Itest)
% $$$   figure
% $$$   spy(Itrain)
% $$$   figure
% $$$   spy(~isnan(Y));
% $$$   figure
% $$$   spy(Itrain&Itest)
% $$$   figure
% $$$   spy((Itrain|Itest)==~isnan(Y));
% $$$   return

  testset_string = 'pattern';
  
 otherwise
  error('Unknown test set requested');
end

switch model
 
 case 1
 
  % 
  % Local GP
  %
  
  %disp('Model: Local GP');
  disp('Model: Local GP with two covfuncs per domain');
  
  % Temporal covariance function
  d = abs(1-(1:length(data.time)));
% $$$   covfunc1 = gp_cov_toeplitz(gp_cov_pp(d,1));
  covfunc1 = gp_cov_toeplitz(gp_cov_sum(gp_cov_pp(d,1), ...
                                        gp_cov_scale(gp_cov_pp(d,1))));
% $$$   theta_temporal = [7];   % length scale
  theta_temporal = [7;   % length scale 1
                    1.0; % magnitude 2
                    3];  % length scale 2
 
  % Spatial covariance function
  [LON,LAT] = meshgrid(data.longitude,data.latitude);
  X = geographic_to_euclidean([LON(:)';LAT(:)']);

  % Use block-Toeplitz structure for the covariance function
  [lat,lon0] = meshgrid(data.latitude,data.longitude(1));
  X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
  D = sqrt(sq_dist(X0,X));
% $$$   covfunc2 = gp_cov_toeplitz_block(gp_cov_pp(D,3));
  covfunc2 = gp_cov_toeplitz_block(gp_cov_sum(gp_cov_pp(D,3), ...
                                              gp_cov_scale(gp_cov_pp(D,3))));

  % Select sea areas
  covfunc2 = gp_cov_select(covfunc2, sea);

% $$$   theta_spatial = [3000];  % length scale
  theta_spatial = [3000;  % length scale 1
                   1.0;   % magnitude 2
                   2000]; % length scale 2
 

  % Initial guess for covariance parameters
  theta_init = [0.5; ...               % total magnitude
                theta_temporal(:); ... % temporal parameters
                theta_spatial(:); ...  % spatial parameters
                0.5]';                 % noise magnitude

  % Prior for the hyperparameters
  a = 1e-3;
  b = 1e-3;
  logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
  dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);
  
  % Noise scaling
  w = 1./sqrt(cosd(LAT));
  noise_scale = repmat(w(sea), [1,size(Y,2)]);

  % Inference
  N_samples = 1000;
  burnin = floor(N_samples/2);
  res = gp_kron(Ytrain, covfunc1, covfunc2, N_samples, theta_init, ...
                logprior_theta, dlogprior_theta, burnin, 'noise_scale', ...
                noise_scale);
  
  F = res.F;
  stdF = sqrt(res.FF - res.F.^2);
  results = res;
  
  %model_string = 'gpkron';
  model_string = 'gpkron-cov2';
  
  
 case 2
  
  %
  % VB PCA
  %
 
  disp('Model: VB PCA');
  
  D = 150;

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
  noise_module = noise_module_isotropic(M, N, 1e-3, 1e-3, ...
                                        'init', 10, ...
                                        'weights', weights);

  % Run VB PCA
  Q = vbfa(D, Ytrain, W_module, X_module, noise_module, ...
           'maxiter', 500, ...
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

  model_string = 'vbpca';
  
 otherwise
  error('Unknown model requested');
end

% Error measures
[LON,LAT] = meshgrid(data.longitude, data.latitude);
weights = cosd(LAT(:));
weights = weights(sea);
weights = repmat(weights, [1, N]);
rmse_train = rmsew(Y(Itrain)-F(Itrain),weights(Itrain));
rmse_test = rmsew(Y(Itest)-F(Itest),weights(Itest));
rmse_zero = rmsew(Y(Itest),weights(Itest));
fprintf('Training WRMSE=%.4f and testing WRMSE=%.4f\n', rmse_train, rmse_test);
fprintf('Testing WRMSE=%.4f with zero predictions\n', rmse_zero);


filename = sprintf(['/home/jluttine/matlab/publications/nips2011/' ...
                    'results_nips2011_mohsst5_%s_%s_%s'], model_string, ...
                   testset_string, datestr(now,'yyyymmdd'));

save(filename, 'F', 'stdF', 'results', 'Ytest', 'Ytrain');
fprintf('Saved results to %s\n', filename);

if nargout < 1
  clear F
  clear stdF
  clear results
end

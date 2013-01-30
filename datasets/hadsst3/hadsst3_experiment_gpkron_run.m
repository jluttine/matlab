
function res = hadsst3_experiment_gpkron_run()

%
% Process data
%

data = hadsst3_load_data();

% Form the data matrix
Y = data.observations;
[M,N] = size(Y);
Obs = ~isnan(Y);


%
% Temporal covariance function
%

d = abs(1-(1:length(data.time)));
covfunc1 = gp_cov_toeplitz(gp_cov_pp(d,1));
theta_temporal = [7];   % length scale 1

%
% Spatial covariance function
%

[LON,LAT] = meshgrid(data.longitude,data.latitude);

X = geographic_to_euclidean([LON(:)';LAT(:)']);

% Use block-Toeplitz structure for the covariance function
[lat,lon0] = meshgrid(data.latitude,data.longitude(1));
X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
D = sqrt(sq_dist(X0,X));
covfunc2 = gp_cov_toeplitz_block(gp_cov_pp(D,3));
theta_spatial = [3000];  % length scale 1

%
% Inference
%

N_samples = 1000;
burnin = floor(N_samples/2);

folder = sprintf('/share/bayes/jluttine/results/hadsst3/gpkron');
mkdir(folder);
mkdir([folder, '/samples']);
datestring = datestr(now,30)
filename = sprintf('%s/results_hadsst3_gpkron_%s', ...
                   folder, ...
                   datestring);
filename_samples = sprintf('%s/samples/samples_hadsst3_gpkron_%s', ...
                           folder, ...
                           datestring);

% Prior
a = 1e-3;
b = 1e-3;
logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);

% Initial guess for covariance parameters
theta_init = [0.5; ...               % total magnitude
              theta_temporal(:); ... % temporal parameters
              theta_spatial(:); ...  % spatial parameters
              0.5]';                 % noise magnitude

% Weights for the noise levels using the respective grid size
w = 1./sqrt(cosd(LAT));

% Run sampler
res = gp_kron(Y, ...
              covfunc1, ...
              covfunc2, ...
              N_samples, ...
              theta_init, ...
              logprior_theta, ...
              dlogprior_theta, ...
              burnin, ...
              'noise_scale2', w, ...
              'filename_samples', filename_samples)

% Save results
save(filename, '-struct', 'res');
disp(['Saved results to ', filename]);



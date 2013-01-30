%function test_gp

%
% Generate data
%


t = 1:1:200;
N = length(t);

if 1
switch 2
 case 1
  covfunc_data = gp_cov_sum(gp_cov_rq(sq_dist(t,t)), ...
                            gp_cov_scale(gp_cov_delta(N)));
  K = covfunc_data([10,    % length scale
                    1.0,   % alpha
                    0.1]); % noise magnitude
  L = chol(K, 'lower');
  y = L*randn(length(L),1);
 case 2
  y = sin(2*pi*0.03*t(:)) + 0.2*randn(N,1);
end

%y = y + 1e3;
end

%yy = sqrt(mean(y.^2))

%
% Construct covariance functions
%

% Covariance function of the signal (RQ for residuals and SE for trend/bias)
switch 3
case 1
    covfunc = @(D2) gp_cov_sum(gp_cov_scale(gp_cov_rq(D2)), ...
                               gp_cov_scale(gp_cov_se(D2)));
    theta_f = [0.5  % RQ: magnitude
               1    % RQ: alpha
               3    % RQ: length scale
               1e1  % SE: magnitude
               100]; % SE: length scale
case 2
    covfunc = @(D2) gp_cov_scale(gp_cov_rq(D2));
    theta_f = [0.5  % RQ: magnitude
               1    % RQ: alpha
               3 ];   % RQ: length scale
case 3
    covfunc = @(D2) gp_cov_scale(gp_cov_se(D2));
    theta_f = [1e1  % SE: magnitude
               10]; % SE: length scale
    
end

covfunc_f = gp_cov_toeplitz(covfunc(sq_dist(t(1),t)));

% Covariance function of the noisy observations
covfunc_noise = gp_cov_scale(gp_cov_delta(N));
theta_noise = [0.4]; % noise magnitude

%
% Perform GP inference
%

% Learn the hyperparameters
covfunc_y = gp_cov_sum(covfunc_f, covfunc_noise);
theta_y = [theta_f; theta_noise];
theta_y = gp_learn(y, covfunc_y, theta_y, ...
                   'maxiter', 1000, ...
                   'checkgrad', true)

% Select the hyperparameters for the signal
n_theta_f = covfunc_f();
theta_f = theta_y(1:n_theta_f);

% Covariance functions for the predictions
th = (min(t)-10):(max(t)+150);
Nh = length(th);
% $$$ D2_yh = sq_dist(t, th);
% $$$ d2_h = zeros(size(th(:)));
% $$$ covfunc_yh = covfunc(sqrt(D2_yh));
% $$$ covfunc_h = covfunc(sqrt(d2_h));
covfunc_yh = covfunc(sq_dist(t, th));
covfunc_h = covfunc(zeros(Nh,1));

% Compute the distribution of the function values
K_y = covfunc_y(theta_y);
K_yh = covfunc_yh(theta_f);
k_h = covfunc_h(theta_f);
[f_mean, f_var] = gp_predict(y, K_y, K_yh, k_h);

% Alex: Show the posterior covariance function

%
% Show results
%

figure
clf
%errorplot(th, f_mean, 2*sqrt(f_var), 'Color', 'r');
errorplot(th, f_mean, 2*(sqrt(f_var)+theta_y(end)), 'Color', 'r');
hold on
plot(t,y,'k+')

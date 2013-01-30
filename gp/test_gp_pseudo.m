function test_gp_pseudo

%
% Generate data
%

t = 1:1:1000;
N = length(t);
y = sin(2*pi/N*t(:));
y = y + 0.1*randn(size(y));
% $$$ t = 1:2:100;
% $$$ y = sin(2*pi*0.03*t(:));
% $$$ y = y + 0.2*randn(size(y));
% $$$ N = length(t);

%
% Construct covariance functions
%

% Distance matrices for the covariance functions
tp = linspace(1,N,20);
D2_pp = sq_dist(tp, tp);
D2_pf = sq_dist(tp, t);
d2_f = zeros(size(t(:)));

% Covariance function of the signal
covfunc = @(D2) gp_cov_scale(gp_cov_se(D2));
covfunc_pseudo = gp_cov_pseudo(gp_cov_scale(gp_cov_jitter(gp_cov_se(D2_pp),1e-6)), ...
                               covfunc(D2_pf), ...
                               covfunc(d2_f));
% $$$ covfunc_pseudo = gp_cov_pseudo(covfunc, ...
% $$$                                {D2_pp}, ...
% $$$                                {D2_fp}, ...
% $$$                                {d2_f});
theta_f = [2  % signal magnitude
           5]; % length-scale

% Covariance function of the noise
covfunc_noise = gp_cov_scale(gp_cov_delta(N));
theta_noise = [2]; % noise magnitude

%
% Perform GP inference
%

% Learn the hyperparameters
[theta_f, theta_noise] = gp_learn_pseudo(y, covfunc_pseudo, theta_f, covfunc_noise, ...
                                         theta_noise, ...
                                         'maxiter', 50, ...
                                         'checkgrad', false);
%theta_f = gp_learn_pseudo(y, covfunc_f, theta_f, U);
theta = [theta_f; theta_noise]

% Covariance functions for the predictions
th = -100:(N+150);
D2_ph = sq_dist(tp, th);
d2_h = zeros(size(th(:)));
covfunc_ph = covfunc(D2_ph);
covfunc_h = covfunc(d2_h);

% Compute the distribution of the function values
K_noise = covfunc_noise(theta_noise);
[K_pp, K_fp, ~] = covfunc_pseudo(theta_f);
K_ph = covfunc_ph(theta_f);
k_h = covfunc_h(theta_f);

% $$$ figure
% $$$ imagesc(K_pp)
% $$$ return

[f_mean, f_var] = gp_predict_pseudo(y, K_noise, K_fp, K_pp, K_ph, k_h);

%
% Show results
%

clf
%errorplot(th, f_mean, 2*sqrt(f_var), 'Color', 'r');
errorplot(th, f_mean, 2*(sqrt(f_var)+theta_noise(end)), 'Color', 'r');
hold on
plot(t,y,'k+')
plot(tp,zeros(size(tp)), 'bo');

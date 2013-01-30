
% [covmatrix] = covfunc(theta)
% [covmatrix, theta] = rand_theta(theta, covmatrix, Y, covfunc,
%                                 logposterior)
% [rand_y] = get_randfunc_y(covmatrix)
% Y = rand_y(covmatrix, Y, I)
% [y, dy] = logprior_theta(theta)
% [y, dy] = loglikelihood_theta(Y, covmatrix)

% COVMATRIX:
% L1/LD1
% L2/LD2
% linsolve1
% linsolve2
% K1
% K2
% logdet1
% logdet2

function res = test_gp_kron(N1,N2,N_samples,seed)

%
% DATA
%

%randn('state', 10)
%rand('state', 10)

if nargin == 0
  N1 = 200;
  N2 = 300;
end
if nargin < 3
  N_samples = 100;
end
if nargin >= 4
  randn('state', seed);
  rand('state', seed);
else
  randn('state', 10);
  rand('state', 10);
end
  

% Inputs
x1 = 1:N1;
x2 = 1:N2;

% Covariance models

d1 = sqrt(sq_dist(x1(1), x1));
covfunc1 = gp_cov_pp(d1,1);
covfunc1 = gp_cov_toeplitz(covfunc1);
theta1 = 10;
K1 = covfunc1(theta1);

d2 = sqrt(sq_dist(x2(1), x2));
covfunc2 = gp_cov_pp(d2,1);
covfunc2 = gp_cov_toeplitz(covfunc2);
theta2 = 10;
K2 = covfunc2(theta2);

% $$$ nz1 = nnz(K1)
% $$$ nz2 = nnz(K2)
% $$$ K = kron(K1,K2);
% $$$ whos
% $$$ t = cputime();
% $$$   LD = ldlchol(K);
% $$$   ldlsolve(LD,randn(N1*N2,1));
% $$$ time = cputime() - t
% $$$ return

% $$$ % Plot gradients
% $$$ [K1, dK1] = covfunc1(theta1);
% $$$ figure
% $$$ subplot(3,1,1)
% $$$ imagesc(K1)
% $$$ subplot(3,1,2)
% $$$ imagesc(dK1{1})
% $$$ subplot(3,1,3)
% $$$ imagesc(dK1{2})
% $$$ return

disp('Generate data..')

% Generate noisy data
a = 1;
s = 1;
Y = a * kronprod(lchol(K1), lchol(K2), randn(N2,N1));
% $$$ Y = a * (lchol(K2) * randn(N2,N1) * lchol(K1)');
save('Y_noiseless', 'Y');
Y = Y + s*randn(N2,N1);
save('Y_noisy', 'Y');

% Missing values
pmv = 0.5;
%Imv = randperm(N2*N1);
%Imv = Imv(1:floor(pmv*(N2*N1)));
Imv = rand(N2,N1) < pmv;
%Imv = (rand(N2,N1) < 0.2);
ind1 = ceil(N1/2) + (1:(2*theta1));
ind2 = ceil(N2/2) + (1:(2*theta2));
Imv(ind2,ind1) = true;
Y(Imv) = nan;

%
% INFERENCE
%

a = 1e-3;
b = 1e-3;
logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);

% Gibbs sampling for Y(missing) and covariance parameters
burnin = floor(N_samples/2);

% Initial guess for covariance parameters
theta_init = [2.0        ... % total magnitude
              2.0*theta1 ... % 1) length scale
              0.5*theta2 ... % 2) length scale
              2.0]';         % noise magnitude

res = gp_kron(Y, covfunc1, covfunc2, N_samples, theta_init, ...
              logprior_theta, dlogprior_theta, burnin)

save(sprintf('/home/jluttine/matlab/gp/results_test_gp_kron_%s', ...
             datestr(now,'yyyymmdd')), ...
     '-struct', 'res');

cl = max(abs(Y(:))) + 0.1;
clim = [-cl cl];

Y(Imv) = cl;

figure()
clf

subplot(2,3,1)
imagesc(Y,clim)

load('Y_noiseless')

subplot(2,3,4)
imagesc(Y,clim)
title('true')

subplot(2,3,5)
imagesc(res.F,clim)
title('mean')

subplot(2,3,6)
F_var = sqrt(res.FF - res.F.*res.F);
clim_var = [0, max(F_var(:))];
imagesc(F_var,clim_var)
title('std')

map_colormap();
cm = colormap();
cm(end,:) = [0.5 0.5 0.5];
colormap(cm);

rmse_F = rmse(Y,res.F)

subplot(2,3,2)
semilogy(res.theta');

subplot(2,3,3)
%lag = 200;
plot(acorr(res.theta(:,burnin:end)'));

figure()
plot_scatterhist(res.theta(:,burnin:end)');


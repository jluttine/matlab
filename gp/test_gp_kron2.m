
% Test a kronecker product covariance matrix GP.  The covariance matrix is a
% sum of three Kronecker-product matrices plus noise.

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

function res = test_gp_kron2(N1,N2,N_samples,seed)

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
  N_samples = 4;
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

covfuncs = cell(2,3);
thetas = cell(2,3);
K = cell(2,3);

% Distances on the two domains
d1 = sqrt(sq_dist(x1(1), x1));
d2 = sqrt(sq_dist(x2(1), x2));

% "Spatial" covariance functions
covfunc = gp_cov_pp(d1,1);
covfunc = gp_cov_toeplitz(covfunc);
covfuncs(1,:) = {covfunc};
thetas(1,:) = {5; 20; 50};

% "Temporal" covariance functions
covfunc = gp_cov_pp(d2,1);
covfunc = gp_cov_toeplitz(covfunc);
covfuncs(2,:) = {covfunc};
thetas(2,:) = {5; 20; 50};

% Compute covariance matrices
for i=1:2
  for j=1:3
    K{i,j} = feval(covfuncs{i,j}, thetas{i,j});
  end
end

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
a(1) = 1;
Y = a(1) * kronprod(lchol(K{1,1}), lchol(K{2,1}), randn(N2,N1));
a(2) = 1;
Y = Y + a(2) * kronprod(lchol(K{1,2}), lchol(K{2,2}), randn(N2,N1));
a(3) = 1;
Y = Y + a(3) * kronprod(lchol(K{1,3}), lchol(K{2,3}), randn(N2,N1));
% $$$ Y = a * (lchol(K2) * randn(N2,N1) * lchol(K1)');
save('Y_noiseless', 'Y');

s = 1;
Y = Y + s*randn(N2,N1);
save('Y_noisy', 'Y');

% Missing values
pmv = 0.5;
%Imv = randperm(N2*N1);
%Imv = Imv(1:floor(pmv*(N2*N1)));
Imv = rand(N2,N1) < pmv;
%Imv = (rand(N2,N1) < 0.2);
% $$$ ind1 = ceil(N1/2) + (1:20);
% $$$ ind2 = ceil(N2/2) + (1:20);
% $$$ Imv(ind2,ind1) = true;
Y(Imv) = nan;

%
% INFERENCE
%

% Prior
%
% If you give non-informative priors for the length scales, then if a
% signal is pruned out the length scales may get arbitrary values and
% Cholesky decompositions fail.  Thus, give some reasonable priors for
% the length scales! In addition, the magnitude parameter can go
% arbitrarily close to zero and prior should not give INF then, or
% preferably, the prior should not be improper.
a = 1;
b = [1  ... % signal 1 magnitude
     0.1 ... % signal 1 spatial length scale
     0.1 ... % signal 1 temporal length scale
     1  ... % signal 2 magnitude
     0.1 ... % signal 2 spatial length scale
     0.1 ... % signal 2 temporal length scale
     1  ... % signal 3 magnitude
     0.1 ... % signal 3 spatial length scale
     0.1 ... % signal 3 temporal length scale
     1]';   % noise magnitude
% $$$ a = 1e-3;
% $$$ b = 1e-3;
logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);

burnin = floor(N_samples/2);

% Initial guess for covariance parameters
theta_init = [2.0             ... % signal 1 magnitude
              2.0*thetas{1,1} ... % signal 1 spatial length scale
              0.5*thetas{2,1} ... % signal 1 temporal length scale
              2.0             ... % signal 2 magnitude
              2.0*thetas{1,2} ... % signal 2 spatial length scale
              0.5*thetas{2,2} ... % signal 2 temporal length scale
              2.0             ... % signal 3 magnitude
              2.0*thetas{1,3} ... % signal 3 spatial length scale
              0.5*thetas{2,3} ... % signal 3 temporal length scale
              2.0]';              % noise magnitude

res = gp_kron(Y, ...
              covfuncs(1,:), ...
              covfuncs(2,:), ...
              N_samples, ...
              theta_init, ...
              logprior_theta, ...
              dlogprior_theta, ...
              burnin, ...
              'hyperparameters', 'ss')

save(sprintf('/home/jluttine/matlab/gp/results_test_gp_kron_%s', ...
             datestr(now,'yyyymmdd')), ...
     '-struct', 'res');

%
% PLOT RESULTS
%

cl = max(abs(Y(:))) + 0.1;
clim = [-cl cl];

Y(Imv) = cl;

figure()
clf

subplot(3,3,1)
imagesc(Y,clim)
title('observations')

load('Y_noiseless')

subplot(3,3,2)
imagesc(Y,clim)
title('true')

subplot(3,3,3)
imagesc(sum(res.F,3),clim)
title('mean sum')

subplot(3,3,4)
imagesc(res.F(:,:,1),clim)
title('mean')

subplot(3,3,5)
imagesc(res.F(:,:,2),clim)
title('mean')

subplot(3,3,6)
imagesc(res.F(:,:,3),clim)
title('mean')

F_var = sqrt(res.FF - res.F.*res.F);
clim_var = [0, max(F_var(:))];

subplot(3,3,7)
imagesc(F_var(:,:,1),clim_var)
title('std')

subplot(3,3,8)
imagesc(F_var(:,:,2),clim_var)
title('std')

subplot(3,3,9)
imagesc(F_var(:,:,3),clim_var)
title('std')

map_colormap();
cm = colormap();
cm(end,:) = [0.5 0.5 0.5];
colormap(cm);

rmse_F = rmse(Y,sum(res.F,3))


figure()

subplot(1,2,1)
semilogy(res.theta');

subplot(1,2,2)
plot(acorr(res.theta(:,burnin:end)'));

figure()
plot_scatterhist(res.theta(:,burnin:end)');


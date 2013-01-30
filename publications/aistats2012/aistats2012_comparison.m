function aistats2012_comparison(D, s, sample_hyperparameters, kronecker, fid)

if nargin < 5 || ~isempty(fid)
  fid = 1;
  fprintf(fid, 'Data size: 10^%g\n', D);
  fprintf(fid, 'Noise level: %g\n', s);
  if kronecker
    fprintf(fid, 'Kronecker approach\n');
  else
    fprintf(fid, 'Standard Cholesky approach\n');
  end
end
N = 10^(D);
N1 = ceil(sqrt(N));
N2 = N1;


if kronecker
%  N_samples = 2;
  N_samples = 1000;
  burnin = min(200,ceil(N_samples/2));
else
%  N_samples = 2;
  N_samples = 100;
  burnin = min(20,N_samples/2);
end

seed = 10;
randn('state', seed);
rand('state', seed);

%
% Covariance models
%

% Inputs
x1 = 1:N1;
x2 = 1:N2;

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

disp('Generate data..')

%
% Generate data
%

a = 1;
%s = 0.3;
Y = a * kronprod(lchol(K1), lchol(K2), randn(N2,N1));
F_true = Y;
save('Y_noiseless', 'Y');
Y = Y + s*randn(N2,N1);
save('Y_noisy', 'Y');
Imv = false(size(Y));
%Imv = [];

folder = '/home/jluttine/matlab/publications/aistats2012/results';
filename = sprintf('%s/results_aistats2012_comparison_D=%g_kron=%d_s=%g_date=%s', ...
                   folder, D, kronecker, s, datestr(now, 'yyyymmdd'));



%
% Prior and initialization
%

alpha = 1e-3;
beta = 1e-3;
logprior_theta = @(theta) sum(gamma_logpdf(theta, alpha, beta));
dlogprior_theta = @(theta) gamma_dlogpdf(theta, alpha, beta);

% Initial guess for covariance parameters
theta_init = [1.0*a      ... % total magnitude
              1.0*theta1 ... % 1) length scale
              1.0*theta2 ... % 2) length scale
              1.0*s]';       % noise magnitude


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Just plot some reconstructions (COMMENT OUT)
%
if true
  % Standard approach
  covfunc_f = gp_cov_kron(gp_cov_scale(covfunc1), covfunc2);
  covfunc_y = gp_cov_sum(covfunc_f, gp_cov_scale(gp_cov_delta(N1*N2)));
  K_f = covfunc_f(theta_init(1:(end-1)));
  K_y = covfunc_y(theta_init);
  k_f = diag(K_f);
  % Compute the variances
  f_var = zeros(size(k_f));
  [L,p,q] = lchol(K_y);
  n2 = 10;%floor(N2/2);
  for n1=1:N1
    n = n1 + (n2-1)*N1;
    z(q) = L\K_f(q,n);
    f_var(n) = k_f(n) - z(:)'*z(:);
  end
  F_std2 = 2*sqrt(reshape(f_var, [N1 N2]));
% $$$   for n=1:size(K_y,2)
% $$$     z(q) = L\K_f(q,n);
% $$$     f_var(n) = k_f(n) - z(:)'*z(:);
% $$$   end
  y = Y(:);
  F_mean = reshape(K_f(:,q)*(L'\(L\y(q))), [N1 N2]);
  res1 = gp_kron(Y, covfunc1, covfunc2, 1, theta_init, logprior_theta, ...
                 dlogprior_theta, 0, 'hyperparameters', 'fixed')
  res10 = gp_kron(Y, covfunc1, covfunc2, 10, theta_init, logprior_theta, ...
                  dlogprior_theta, 0, 'hyperparameters', 'fixed')
% $$$   res100 = gp_kron(Y, covfunc1, covfunc2, 100, theta_init, logprior_theta, ...
% $$$                    dlogprior_theta, 0, 'hyperparameters', 'fixed')
  res10.std2 = 2 * sqrt(res10.FF - res10.F.^2);
% $$$   res100.std2 = 2 * sqrt(res100.FF - res100.F.^2);
  figure
  style = {'ErrorLineStyle', ':', ...
           'ErrorFillAlpha', 0};
  plot(F_true(:,n2), 'k');
  hold on
  plot(Y(:,n2), '+k');
  plot(F_mean(:,n2), 'r-.');
% $$$   errorplot(res100.F(:,n2), res100.std2(:,n2), 'Color', 'b', 'ErrorEdgeColor', ...
% $$$             'b', style{:});
  plot(res10.F(:,n2), 'b--');
% $$$   plot(res1.F(:,n2), 'g');
  errorplot(F_mean(:,n2), F_std2(:,n2), 'Color', 'r', 'LineStyle', '-.', ...
            'ErrorEdgeColor', 'r', style{:});
  errorplot(res10.F(:,n2), res10.std2(:,n2), 'Color', 'b', 'LineStyle', ...
            '--', 'ErrorEdgeColor', 'b', style{:});
  set(gca, 'xlim', [1 N1]);
  hl = legend({'True function', 'Observations', 'Exact posterior', ...
               '10 samples'}, 'Location', 'NorthWest');
  %set(gca, 'Box', 'off' );
  set(gca, 'FontSize', 10);
  set(hl, 'FontSize', 8);
  hx = xlabel('Time (t)');
  set(hx, 'FontSize', 10);
  set(gca, 'Box', 'off');
  set(gca, 'xlim', [10 25])
  set_figure_size(13,11);
  %set_subplot_positions(gca, 1, 1, [0.15 0.1 0.1 0.15], [0 0]);
  print('-depsc2', '/home/jluttine/papers/gpkron/figures/fig_artificial_comparison');
  return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if kronecker
  % Kronecker trick
  if sample_hyperparameters
    time = cputime();
    res = gp_kron(Y, covfunc1, covfunc2, N_samples, theta_init, ...
                  logprior_theta, dlogprior_theta, burnin)
    time = cputime() - time;
    theta = res.theta;
    fprintf(fid, 'Time for inference of theta: %g seconds\n', time);
  else
    time = cputime();
    res = gp_kron(Y, covfunc1, covfunc2, 10, theta_init, logprior_theta, ...
                  dlogprior_theta, burnin, 'hyperparameters', 'fixed')
    time = (cputime() - time) / (10+1)
    fprintf(fid, 'Time for inference of F: %g seconds\n', time);
  end
else
  % Standard approach
  covfunc_f = gp_cov_kron(gp_cov_scale(covfunc1), covfunc2);
  covfunc_y = gp_cov_sum(covfunc_f, gp_cov_scale(gp_cov_delta(N1*N2)));
  
  if sample_hyperparameters
    time = cputime();
    theta = gp_mcmc(Y(:), covfunc_y, N_samples, theta_init, logprior_theta, ...
                    dlogprior_theta);
    fprintf(fid, 'Time for inference of theta: %g seconds\n', cputime()-time);
    time = cputime() - time;
  else
    K_f = covfunc_f(theta_init(1:(end-1)));
    K_y = covfunc_y(theta_init);
    k_f = diag(K_f);
    time = cputime();
    disp('Perform inference')
    % Compute only the variances for simplicity
    f_var = zeros(size(k_f));
    %invK_y = spinv(K_y);
    [L,p,q] = lchol(K_y);
    for n=1:size(K_y,2)
      z(q) = L\K_f(q,n);
      f_var(n) = k_f(n) - z(:)'*z(:);
      %f_var(n) = k_f(n) - K_f(:,n)'*(invK_y*K_f(:,n));
    end
    %[f_mean, f_var] = gp_predict(Y(:), K_y, K_f, k_f);
    fprintf(fid, 'Time for inference of F: %g seconds\n', cputime()-time);
    time = cputime() - time
    % [LD,p,q] = ldlchol(K);
  end
    
end

if sample_hyperparameters
  save(filename, 'theta', 'time');

  % Effective number of samples:

  ess = zeros(size(theta,1),1);
  for m=1:size(theta,1)
    ess(m) = mcmc_ess_acorr(theta(m,burnin:end));
  end
  ess
  ess = min(ess);
  fprintf(fid, 'Number of samples after burn-in: %d\n', N_samples-burnin+1);
  fprintf(fid, 'CPU time per sample: %.3f seconds\n', time/N_samples);
  fprintf(fid, 'Effective number of samples: %.2f\n', ess);
  time_ess = time * (N_samples-burnin+1)/N_samples;
  fprintf(fid, 'CPU time per effective sample: %.3f seconds\n', time_ess/ess);

% $$$   figure()
% $$$   clf
% $$$ 
% $$$   cl = max(abs(Y(:)));
% $$$   clim = [-cl cl];
% $$$ 
% $$$   subplot(2,2,1)
% $$$   imagesc(Y,clim)
% $$$ 
% $$$   load('Y_noiseless')
% $$$ 
% $$$   subplot(2,2,2)
% $$$   imagesc(Y,clim)
% $$$   title('true')
% $$$ 
% $$$   subplot(2,2,3)
% $$$   semilogy(theta');
% $$$ 
% $$$   subplot(2,2,4)
% $$$   %lag = 200;
% $$$   plot(acorr(theta(:,burnin:end)'));
% $$$ 
% $$$ % $$$ subplot(2,3,5)
% $$$ % $$$ res.F
% $$$ % $$$ imagesc(res.F,clim)
% $$$ % $$$ title('mean')
% $$$ % $$$ 
% $$$ % $$$ subplot(2,3,6)
% $$$ % $$$ F_var = sqrt(res.FF - res.F.*res.F);
% $$$ % $$$ clim_var = [0, max(F_var(:))];
% $$$ % $$$ imagesc(F_var,clim_var)
% $$$ % $$$ title('std')
% $$$ 
% $$$   map_colormap();
% $$$ % $$$ cm = colormap();
% $$$ % $$$ cm(end,:) = [0.5 0.5 0.5];
% $$$ % $$$ colormap(cm);
% $$$ 
% $$$ % $$$ rmse_F = rmse(Y,res.F)
% $$$ 
% $$$   figure()
% $$$   plot_scatterhist(theta(:,burnin:end)');
end

fprintf(fid, '\n');


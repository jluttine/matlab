function  theta = gp_mcmc(y, covfunc, N_samples, theta_init, logprior_theta, ...
                          dlogprior_theta)

%
% Posterior distribution
%

get_logpdf = @(theta) get_logpdf_theta(theta);
% $$$ get_logpdf = @(theta) (@() (gp_loglikelihood(y,covfunc,theta,'sparse2dense',0) ...
% $$$                             + logprior_theta(theta)));
get_dlogpdf = @(theta) (@() nan());
% $$$ get_logpdf = @(theta) (@() (gp_loglikelihood(y, covfunc, theta) + ...
% $$$                             logprior_theta(theta)));
% $$$ get_dlogpdf = @(theta) (@() nan());

%
% Transform to log-scale
%

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
get_logpdf = @(logtheta) get_logpdf_transformed(exp(logtheta), ...
                                               get_logpdf, ...
                                               sum(logtheta));
get_dlogpdf = @(logtheta) get_dlogpdf_transformed(exp(logtheta), ...
                                                  get_dlogpdf, ...
                                                  diag(exp(logtheta)), ...
                                                  ones(size(logtheta)));

logtheta_init = log(theta_init);

%
% Run slice sampler
%

theta = nan(length(theta_init), N_samples);
rand_theta = mcmc_init_slicesampling(logtheta_init, get_logpdf);
for n=1:N_samples
  t = cputime();
  theta(:,n) = exp(rand_theta());
  fprintf('Theta sample %d (%.3f seconds)\n', n, cputime()-t);
end

% $$$   function f_theta = func_theta(theta)
% $$$   % Pre-compute stuff.
% $$$   f_theta.loglikelihood = gp_loglikelihood(y, covfunc, theta);
% $$$   f_theta.theta = theta;
% $$$ % $$$   K = covfunc(theta);
% $$$ % $$$   LD = ldlchol(K);
% $$$ % $$$   f_theta.y_invK_y = y'*linsolve_ldlchol(LD,Y);
% $$$ % $$$   f_theta.logdetK = logdet_ldlchol(LD);
% $$$ % $$$   f_theta.theta = theta;
% $$$   end

  function func = get_logpdf_theta(theta)
  logpdf = gp_loglikelihood(y,covfunc,theta,'sparse2dense',0) + ...
           logprior_theta(theta);
  func = @() logpdf;
  end
end



% $$$ function [f_theta, df_theta] = func_theta_transformed(theta_transformed, ...
% $$$                                                   theta, func_theta, varargin)
% $$$ if nargout <= 1
% $$$   f_theta = func_theta(theta, varargin{:});
% $$$   f_theta.theta_transformed = theta_transformed;
% $$$ else
% $$$   [f_theta, df_theta] = func_theta(theta, varargin{:});
% $$$   f_theta.theta_transformed = theta_transformed;
% $$$   df_theta.theta_transformed = theta_transformed;
% $$$ end
% $$$ end

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


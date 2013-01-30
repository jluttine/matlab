% GP_LEARN_PSEUDO - Learns the hyperparameters of a Gaussian process
%                   using pseudo inputs.
%
% Usage:
%
%   [THETA_F, THETA_NOISE, LOGLIKELIHOOD] = ...
%     GP_LEARN_PSEUDO(Y, COVFUNC_PSEUDO, THETA0_F, COVFUNC_NOISE, THETA0_NOISE)
%
% Y              : Nx1 vector of observations
% COVFUNC_PSEUDO : A special covariance function, see GP_COV_PSEUDO
% THETA0_F       : Initial parameter values for COVFUNC_PSEUDO
% COVFUNC_NOISE  : A diagonal noise covariance function
% THETA0_NOISE   : Initial parameter values for COVFUNC_NOISE
%
% THETA_F       : Optimized parameter values for COVFUNC_PSEUDO
% THETA_NOISE   : Optimized parameter values for COVFUNC_NOISE
% LOGLIKELIHOOD : The loglikelihood lower bound at the optimum
%
% Optional parameters can be given as
%
%   [...] = GP_LEARN_PSEUDO(..., 'PARAMETER', VALUE)
%
% where the possible parameters are
%
% 'MAXITER' : Maximum number of iterations (default: 100)
% 'CHECKGRAD' : Numerical check of the gradients (default: false)
%
% See also GP_PREDICT_PSEUDO, GP_LEARN, GP_COV_PSEUDO,
% GP_LOGLIKELIHOOD_PSEUDO.

% Last modified 2010-01-27
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function [theta_f, theta_noise, loglikelihood] = gp_learn_pseudo(y, ...
                                                  covfunc_pseudo, ...
                                                  theta_f, ...
                                                  covfunc_noise, ...
                                                  theta_noise, ...
                                                  varargin)

options = struct('maxiter', 100, ...
                 'checkgrad', false);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

n_theta_f = numel(theta_f);
theta = [theta_f(:); theta_noise(:)];

if options.checkgrad
  mycheckgrad(@cost, log(theta), 1e-6);
end


[logtheta, logpdfs] = minimize(log(theta), @cost, options.maxiter);
loglikelihood = -min(logpdfs);
theta = exp(logtheta);

  function [c, dc] = cost(logtheta)
  theta = exp(logtheta);
  theta_f = theta(1:n_theta_f);
  theta_noise = theta((n_theta_f+1):end);
  [c, dc_f, dc_noise] = gp_loglikelihood_pseudo(y, covfunc_pseudo, theta_f, ...
                                                covfunc_noise, theta_noise);
  c = -c;
  dc = -[dc_f; dc_noise] .* theta;
  end

end

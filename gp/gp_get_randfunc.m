function func = gp_get_randfunc(covfunc, theta, )

% gp_

% Takes:
% - a function handle covfunc which returns a covmat-struct (its
%   specification is up to the user)
% - a function handle get_randfunc_y which returns sampler function
%   handle rand_y
% - a function handle rand_theta?
% - function handles loglikelihood and logprior?
%
% [covmatrix] = covfunc(theta)
% [covmatrix, theta] = rand_theta(theta, covmatrix, Y, covfunc,
%                                 loglikelihood, logprior)
% [rand_y] = get_randfunc_y(covmatrix)
% Y = rand_y(covmatrix, Y, I)
% [y, dy] = logprior(theta)
% [y, dy] = loglikelihood(Y, covmatrix)

% this system tries to fulfill the requirements of the Gibbs sampler
% for GP + FA model - don't make it more general.. the main purpose is
% to model the residuals and sample missing observations

% Options:
%
% - proposal distribution for theta

% Initialize stuff with the given theta initialization
covmatrix = covfunc(theta);

func = @randfunc;

% returns y-sampler and some solver for FA model?
% covmatrix_linsolver = [];

  function Y = randfunc(Y,I)
  % Reconstruct Y
  Y = rand_y(covmatrix, Y, I);
  % Sample theta
  % logposterior = @(covmatrix, theta) (loglikelihood(covmatrix,Y) + ...
  %                                    logprior(theta));
  [covmatrix,theta_new] = rand_theta(theta, covmatrix, Y, covfunc, ...
                                     loglikelihood, logprior);
  if ~isequal(theta,theta_new)
    rand_y = get_randfunc_y(covmatrix);
  end
  end
  
end
  
  

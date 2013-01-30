% NOISE_MODULE_ISOTROPIC  -  Isotropic noise module for variational
%                            Bayesian factor models
%
% MODULE = NOISE_MODULE_ISOTROPIC(M, N, ...)
%
% M : Rows in the data matrix
% N : Columns in the data matrix
%
% Optional parameters:
% 'prior' : A struct containing parameters for the prior Gamma
%           distribution of tau:
%             'a_tau' : Shape parameter (default: 1E-5)
%             'b_tau' : Scale parameter (default: 1E-5)
% 'init'  : A struct containing initial values for the posterior
%           Gamma distribution of tau:
%             'a_tau' : Shape parameter (default: prior)
%             'b_tau' : Scale parameter (default: prior)

% Last modified 2011-08-16
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function module = noise_module_isotropic(M, N, varargin)

%
% Parse options
%
options = struct('init', struct(), ...
                 'prior', struct(), ...
                 'update_tau', 1);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

%
% Parse prior
%
prior = struct('a_tau', 1e-5, ...
               'b_tau', 1e-5);
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

%
% Initialization (default from the prior)
%
init = prior;
[init, errmsg] = argparse(init, options.init);
error(errmsg);
a_tau = NaN;
b_tau = NaN;
a_tau(:) = init.a_tau;
b_tau(:) = init.b_tau;
tau = a_tau / b_tau;
logtau = psi(a_tau) - log(b_tau);

module.get_struct = @get_struct;
module.update = @update;

  function initialize()
  error('This is for bacward compatibility only!')
  end

  function S = get_struct()
  S.module = 'noise_module_isotropic';
  S.Tau = repmat(tau, M, N);
  S.posterior = struct('a_tau', a_tau, ...
                       'b_tau', b_tau);
  S.prior = prior;
  S.init = init;
  S.options = options;
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  % OBS is the number of observations for each element of E2.
  
  % <(Y-WX)^2>
  e2 = sum(E2(:));
%  e2 = E2(:)' * weights(:);
% $$$   e2 = sum(E2(Obs));
  
  %
  % Update tau
  %
  
  MN_obs = sum(Obs(:));
  
  if index_selected(iter, options.update_tau)
    
    a_tau = prior.a_tau + 0.5*MN_obs;
    b_tau = prior.b_tau + 0.5*e2;
  
    tau = a_tau / b_tau;
    logtau = psi(a_tau) - log(b_tau);
    
    %
    % Compute Kullback-Leibler divergence KL(q(tau)||p(tau))
    %
    
    % <log p(Tau)>
    logpdf_p = gamma_logpdf(tau, prior.a_tau, prior.b_tau, logtau);
    
    % <log q(Tau)>
    logpdf_q = -gamma_entropy(a_tau,b_tau);
    
    % KL(q||p)
    KL_Tau = logpdf_q - logpdf_p;
    
  else
    logtau = log(tau);
    KL_Tau = 0;
  end
  
  % Compute expectations of weighted Tau
  Tau = repmat(tau, M, N);
  LogTau = repmat(logtau, M, N);
% $$$   Tau = tau * weights;
% $$$   LogTau = logtau + log(weights);
  
  end

end

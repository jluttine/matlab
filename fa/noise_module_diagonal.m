% NOISE_MODULE_DIAGONAL  -  Diagonal noise module for variational
%                           Bayesian factor models
%
% MODULE = NOISE_MODULE_DIAGONAL(M, N, ...)
%
% M : Rows in the data matrix
% N : Columns in the data matrix
%
% Optional parameters:
% 'separation' : Tells whether the noise level varies between rows or
%                columns of the data matrix.
%                  'columns' : Different noise level for each column; 
%                              1xN noise vector TAU
%                  'rows'    : Different noise level for each row;
%                              Mx1 noise vector TAU
% 'prior' : A struct containing parameters for the prior Gamma
%           distributions of TAU(J):
%             'a_tau' : Shape parameter; Mx1 or 1xN vector
%                       (default: 1E-5) 
%             'b_tau' : Scale parameter; Mx1 or 1xN vector 
%                       (default: 1E-5)
% 'init'  : A struct containing initial values for the posterior
%           Gamma distributions of TAU(J):
%             'a_tau' : Shape parameter; Mx1 or 1xN vector 
%                       (default: prior) 
%             'b_tau' : Scale parameter; Mx1 or 1xN vector 
%                       (default: prior) 

% Last modified 2011-08-22
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function module = noise_module_diagonal(M,N,varargin)

% Parse options
options = struct('init', struct(), ...
                 'prior', struct(), ...
                 'separation', 'rows', ...
                 'update_tau', 1);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Parse prior
prior = struct('a_tau', 1e-5, ...
               'b_tau', 1e-5);
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

% Initialize (default from the prior)
init = prior;
[init, errmsg] = argparse(init, options.init);
error(errmsg);
switch options.separation
 case 'rows'
  a_tau = nan(M,1);;
  b_tau = nan(M,1);
 case 'columns'
  a_tau = nan(1,N);
  b_tau = nan(1,N);
 otherwise
  error('Unknown separation');
end
a_tau(:) = init.a_tau(:);
b_tau(:) = init.b_tau(:);
tau = a_tau ./ b_tau;
logtau = psi(a_tau) - log(b_tau);

module.update = @update;
module.get_struct = @get_struct;

  function S = get_struct()
  S.module = 'noise_module_isotropic';
  S.Tau = bsxfun(@plus, tau, zeros(M,N));
  S.posterior = struct('a_tau', a_tau, ...
                       'b_tau', b_tau);
  S.prior = prior;
  S.init = init;
  S.options = options;
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  % OBS is the number of observations for each element of E2.
  
  % <(Y-WX)^2>
  switch options.separation
   case 'rows'
    e2 = sum(E2,2);
    obs = sum(Obs,2);
   case 'columns'
    e2 = sum(E2,1);
    obs = sum(Obs,1);
   otherwise
    error('Unknown separation');
  end
    
  %
  % Update tau
  %
  
  if index_selected(iter, options.update_tau)
    
    a_tau(:) = prior.a_tau(:) + 0.5*obs(:);
    b_tau(:) = prior.b_tau(:) + 0.5*e2(:);
  
    tau = a_tau ./ b_tau;
    logtau = psi(a_tau) - log(b_tau);
    
    %
    % Compute Kullback-Leibler divergence KL(q(tau)||p(tau))
    %
    
    % <log p(Tau)>
    logpdf_p = sum(gamma_logpdf(tau(:), prior.a_tau(:), prior.b_tau(:), logtau(:)));
    
    % <log q(Tau)>
    logpdf_q = -sum(gamma_entropy(a_tau,b_tau));
    
    % KL(q||p)
    KL_Tau = logpdf_q - logpdf_p;
    
  else
    logtau = log(tau);
    KL_Tau = 0;
  end
  
  % Compute expectations of weighted Tau
  O = zeros(size(E2));
  Tau = bsxfun(@plus, tau, O);
  LogTau = bsxfun(@plus, logtau, O);
  
  end

end

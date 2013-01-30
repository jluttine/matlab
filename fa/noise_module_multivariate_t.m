% NOISE_MODULE_MULTIVARIATE_T  -  Multivariate Student-t noise module for
%                                 variational Bayesian factor models
%
% MODULE = NOISE_MODULE_MULTIVARIATE_T(M, N, ...)
%
% This module models only the robustness of each sample but not the noise
% level.  Use in combination with other non-robust noise modules.
%
% M : Rows in the data matrix
% N : Columns in the data matrix
%
% Optional parameters:
% 'samples' : Tells whether the multivariate samples are rows or columns
%             of the data matrix.
%               'columns' : Samples are the Mx1 column vectors
%               'rows'    : Samples are the 1XN row vectors
% 'init'    : A struct containing initial values for the degrees of
%             freedom and the posterior Gamma distributions of U(J):
%               'nu'  : Degrees of freedom (default: 1)
%               'a_u' : Shape parameters, Mx1 or 1xN vector (default: 0.5)
%               'b_u' : Scale parameter, Mx1 or 1xN vector (default: 0.5)

% Last modified 2011-08-22
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function module = noise_module_multivariate_t(M,N,varargin)


options = struct('init', struct(), ...
                 'samples', 'columns', ...
                 'update_nu', true, ...
                 'update_u', true);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

switch options.samples
 case 'columns'
  dim = 1;
 case 'rows'
  dim = 2;
 otherwise
  error('Unknown value for ''samples''');
end

% Allocate memory
nu = NaN;
switch options.samples
 case 'columns'
  a_u = nan(1,N);
  b_u = nan(1,N);
 case 'rows'
  a_u = nan(M,1);
  b_u = nan(M,1);
 otherwise
  error('Unknown value for ''samples''');
end

% Parse custom initialization
init = struct('a_u', 1/2, ...
              'b_u', 1/2, ...
              'nu', 1);
[init, errmsg] = argparse(init, options.init);
error(errmsg);

% Initialize
nu(:) = init.nu;
a_u(:,:) = init.a_u;
b_u(:,:) = init.b_u;
u = a_u./b_u;
logu = psi(a_u) - log(b_u);

module.update = @update;
module.get_struct = @get_struct;

  function S = get_struct()
  S.module = 'noise_module_multivariate_t';
  S.Tau = bsxfun(@plus, u, zeros(M,N));
  S.posterior = struct('nu', nu, ...
                       'a_u', a_u, ...
                       'b_u', b_u);
  S.init = init;
  S.options = options;
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  % OBS is the number of observations for each element of E2.
  
  [M,N] = size(E2);
  
  %
  % Update nu
  %
  
  if index_selected(iter, options.update_nu)
    nu = t_ml(sum(E2,dim), nu, sum(Obs,dim));
  end

  %
  % Update u
  %
  
  if index_selected(iter, options.update_u)
    a_u = (nu + sum(Obs,dim)) / 2;
    b_u = (nu + sum(E2,dim)) / 2;
    u = a_u./b_u;
    logu = psi(a_u) - log(b_u);
  end
  
  %
  % KL-term: <log q(U)> - <log p(U)>
  %
  
  KL_u = - sum(gamma_entropy(a_u,b_u)) ...
         - sum(gamma_logpdf(u,nu/2,nu/2,logu));

  
  Tau = bsxfun(@times, u, ones(M,N));
  LogTau = bsxfun(@plus, logu, zeros(M,N));
  KL_Tau = KL_u;
  
  
  end

end

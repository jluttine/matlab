% module = noise_module_multivariate_t(noise_module, varargin)
%
% 'samples' can be 'columns' or 'rows'

function module = noise_module_multivariate_t(varargin)


options = struct('init', struct(), ...
                 'samples', 'columns', ...
                 'update_nu', true, ...
                 'update_u', true);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

a_u = [];
b_u = [];
u = [];
nu = [];

switch options.samples
 case 'columns'
  dim = 1;
 case 'rows'
  dim = 2;
 otherwise
  error('Unknown value for ''samples''');
end

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

  function [Tau] = initialize(M, N)

  % Allocate memory
  nu = NaN;
  switch options.samples
   case 'columns'
    u = nan(1,N);
   case 'rows'
    u = nan(M,1);
   otherwise
    error('Unknown value for ''samples''');
  end
  a_u = nan(M,N);
  b_u = nan(M,N);

  % Parse custom initialization
  init = struct('u', 1, ...
                'nu', 1);
  [init, errmsg] = argparse(init, options.init);
  error(errmsg);
  
  % Initialize
  u(:) = init.u;
  nu(:) = init.nu;
  Tau = bsxfun(@times, u, ones(M,N));
  end
  
  function S = get_struct()
  % TODO: not ready yet..
  S.module = 'noise_module_multivariate_t';
  S.nu = nu;
  S.u = u;
  S.a_u = a_u;
  S.b_u = b_u;
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

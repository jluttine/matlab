function module = noise_module_diagonal(varargin)

options = struct('init', struct(), ...
                 'prior', struct(), ...
                 'separation', 'rows', ...
                 'update_tau', 1);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

init = struct('tau', 1);
[init, errmsg] = argparse(init, options.init);
error(errmsg);

prior = struct('a_tau', 1e-5, ...
               'b_tau', 1e-5);
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

tau = [];
a_tau = [];
b_tau = [];
weights = [];

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

  function [Tau] = initialize(M,N)
  switch options.separation
   case 'rows'
    tau = nan(M,1);
    a_tau = nan(M,1);;
    b_tau = nan(M,1);
   case 'columns'
    tau = nan(1,N);
    a_tau = nan(1,N);
    b_tau = nan(1,N);
   otherwise
    error('Unknown separation');
  end
    
  tau(:) = init.tau(:);
  Tau = bsxfun(@plus, tau, zeros(M,N));
  end
  
  function S = get_struct()
  % TODO: not ready yet..
  S.module = 'noise_module_isotropic';
  S.tau = tau;
  S.a_tau = a_tau;
  S.b_tau = b_tau;
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

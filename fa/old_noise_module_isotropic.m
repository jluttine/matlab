function module = noise_module_isotropic(varargin)

options = struct('init', struct(), ...
                 'prior', struct(), ...
                 'weights', 1, ...
                 'update_tau', 1);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

%init = struct('tau', 1);
%[init, errmsg, remopts] = argparse(init, options.init);
%error(errmsg);
tau = 1;
if isfield(options.init, 'tau')
  tau(:,:) = options.init.tau;
end

prior = struct('a_tau', 1e-5, ...
               'b_tau', 1e-5);
[prior, errmsg] = argparse(prior, options.prior);
error(errmsg);

tau = [];
a_tau = [];
b_tau = [];
weights = [];

M = nan;
N = nan;

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

  function [Tau] = initialize(M_rows,N_cols)
  M = M_rows;
  N = N_cols;
  tau = NaN;
  a_tau = [];
  b_tau = [];
  weights = bsxfun(@times, ones(M,N), options.weights);
  tau(:) = init.tau;
  Tau = tau*weights;
  end
  
  function S = get_struct()
  % TODO: not ready yet..
  S.module = 'noise_module_isotropic';
  S.Tau = tau.*ones(M,N);
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
  e2 = E2(:)' * weights(:);
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
  Tau = tau * weights;
  LogTau = logtau + log(weights);
  
  end

end

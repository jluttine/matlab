function module = noise_module_fixed(M, N, tau)

%Tau = bsxfun(@times, tau, ones(M,N));

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

  function S = get_struct()
  S.module = 'noise_module_fixed';
  S.tau = tau;
  S.Tau = bsxfun(@times, tau, ones(M,N));
  S.prior = [];
  S.init = [];
  S.options = [];
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  I = ones(M,N);
  Tau = bsxfun(@times, tau, I);
  LogTau = bsxfun(@times, log(tau), I);
  KL_Tau = 0;

  end

end


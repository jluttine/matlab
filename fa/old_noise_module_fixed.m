function module = noise_module_fixed(tau)

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

  function [Tau] = initialize(M,N)
  Tau = bsxfun(@times, tau, ones(M,N));
  end

  function S = get_struct()
  S.module = 'noise_module_fixed';
  S.tau = tau;
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  I = ones(size(Obs));
  Tau = bsxfun(@times, tau, I);
  LogTau = bsxfun(@times, log(tau), I);
  KL_Tau = 0;

  end

end


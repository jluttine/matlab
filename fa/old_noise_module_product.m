function module = noise_module_product(noise1, noise2)

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

Tau1 = [];
Tau2 = [];

  function [Tau] = initialize(M, N)
  Tau1 = noise1.initialize(M,N);
  Tau2 = noise2.initialize(M,N);
  Tau = Tau1 .* Tau2;
  end
  
  function S = get_struct()
  S.module = 'noise_module_product';
  S.noise1 = noise1.get_struct();
  S.noise2 = noise2.get_struct();
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  [Tau1, LogTau1, KL_Tau1] = noise1.update(iter, Tau2 .* E2, Obs);
  [Tau2, LogTau2, KL_Tau2] = noise2.update(iter, Tau1 .* E2, Obs);
  
  Tau = Tau1 .* Tau2;
  LogTau = LogTau1 + LogTau2;
  KL_Tau = KL_Tau1 + KL_Tau2;
  
  end

end


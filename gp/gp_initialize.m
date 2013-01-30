
function gpstruct = gp_initialize(covstruct)

% Initialize regular GP and use different function for kron and sparse-approx?

covstruct;
gpstruct.get_posterior_sampler = @get_posterior_sampler;

  function sampler = get_posterior_sampler(logtheta)
 
  end
  
  function [mu,s2] = evaluate_posterior_mean(y)
  % should multiply with noiseless covstruct.. :/
  mu = covstruct.multiply( covstruct.linsolve(y) );
  if nargout >= 2
    %
  end
  end
  
end

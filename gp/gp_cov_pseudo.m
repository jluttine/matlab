% COVFUNC = GP_COV_PSEUDO(COVFUNC_HANDLE, X_PP, X_PF, X_F)
% COVFUNC = GP_COV_PSEUDO(COVFUNC_PP, COVFUNC_PF, COVFUNC_F)
%
% [K_PP, K_PF, K_F, DK_PP, DK_PF, DK_F] = COVFUNC(THETA)
%
% NOTE! This function is fundamentally different from all other covariance
% functions GP_COV_*. The returned function COVFUNC is used differently than
% in other covariance functions. Thus, one must NOT pass the returned
% COVFUNC to other GP_COV_* covariance functions!

% Last modified 2010-01-21
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_pseudo(covfunc_handle, x_pp, x_pf, x_f)

if nargin == 4

  if iscell(x_pp)
    covfunc_pp = covfunc_handle(x_pp{:});
  else
    covfunc_pp = covfunc_handle(x_pp);
  end
  
  if iscell(x_pf)
    covfunc_pf = covfunc_handle(x_pf{:});
  else
    covfunc_pf = covfunc_handle(x_pf);
  end

  if iscell(x_f)
    covfunc_f = covfunc_handle(x_f{:});
  else
    covfunc_f = covfunc_handle(x_f);
  end
  
elseif nargin == 3
  
  covfunc_pp = covfunc_handle;
  covfunc_pf = x_pp;
  covfunc_f = x_pf;
  
else
  
  error('Wrong number of inputs');
  
end

covfunc = @get_covariance;

  function [K_pp, K_pf, k_f, dK_pp, dK_pf, dk_f] = get_covariance(theta)
  
  if nargin == 0
    error('Not yet implemented');
  end
  
  if nargout >= 4
    [K_pp, dK_pp] = covfunc_pp(theta);
  else
    K_pp = covfunc_pp(theta);
  end
  
  if nargout >= 5
    [K_pf, dK_pf] = covfunc_pf(theta);
  elseif nargout >= 2
    K_pf = covfunc_pf(theta);
  end
  
  if nargout >= 6
    [k_f, dk_f] = covfunc_f(theta);
  elseif nargout >= 3
    k_f = covfunc_f(theta);
  end
  
  end
  
end
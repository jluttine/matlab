function [func, fx] = mcmc_init_metropolishastings(x_init, get_logpdf, q, ...
                                                  logpdf_q, func_x)

if nargin < 5
  func_x = @(x) x;
end

x_current = x_init;
fx_current = func_x(x_current);
logpdf_current = get_logpdf(fx_current);

func = @metropolishastings;
if nargout >= 2
  fx = fx_current;
end

  function [x, fx] = metropolishastings(varargin)
  
  % Draw proposal
  x_proposal = q(x_current);
  fx_proposal = func_x(x_proposal);
  logpdf_proposal = get_logpdf(fx_proposal);
  lp_proposal = logpdf_proposal(varargin{:});

  % Compare to the current value
  lp_current = logpdf_current(varargin{:});
  r = (lp_proposal - logpdf_q(x_proposal, x_current)) - ...
      (lp_current - logpdf_q(x_current, x_proposal));
  if log(rand) < r
    % Accept
    %disp('accept')
    x_current = x_proposal;
    fx_current = fx_proposal;
    logpdf_current = logpdf_proposal;
  else
    % Reject (do nothing to the current values)
    %disp('reject')
  end
  x = x_current;
  fx = fx_current;
  end
end


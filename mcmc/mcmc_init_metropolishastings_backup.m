function func = mcmc_init_metropolishastings(x_init, get_logpdf, q, logpdf_q, ...
                                             func_x)

if nargin < 6
  samplefunc = [];
end
x_current = x_init;
fx_current = func_x(x_current);
logpdf_current = get_logpdf(x_current, fx_current);

func = @metropolishastings;

  function [x, fx] = metropolishastings(varargin)
  
  %figure
  %rmse = sqrt(mean(mean((varargin{1}-varargin{2}).^2)))
  %imagesc(varargin{1})
  %map_colormap();
  %error('jou')
  
  % Draw proposal
  x_proposal = q(x_current);
  fx_proposal = func_x(x_proposal);
  logpdf_proposal = get_logpdf(x_proposal, fx_proposal);
  lp_proposal = logpdf_proposal(varargin{:});

  % Compare to the current value
  lp_current = logpdf_current(varargin{:});
  r = (lp_proposal - logpdf_q(x_proposal, x_current)) - ...
      (lp_current - logpdf_q(x_current, x_proposal));
  if log(rand) < r
    % Accept
    disp('accept')
    fx_current = fx_proposal;
    x_current = x_proposal;
    logpdf_current = logpdf_proposal;
  else
    % Reject (do nothing to the current values)
    disp('reject')
  end
  x = x_current;
  fx = fx_current;
  end
end


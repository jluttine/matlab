% [func, fx] = mcmc_init_hamiltonian(x_init, get_logpdf, get_dlogpdf, ...
% epsilon, L, func_x)
%
% epsilon is the parameter to the exponential distribution for sampling
% the momentum
%
% L is the number of simulation steps
function [func, fx] = mcmc_init_hamiltonian(x_init, get_logpdf, get_dlogpdf, ...
                                            epsilon, L, func_x)

if nargin < 6
  func_x = @func_x_default;
end

x_current = x_init;
[fx_current, dfx_current] = func_x(x_current);
logpdf_current = get_logpdf(fx_current);
dlogpdf_current = get_dlogpdf(dfx_current);

func = @hamiltonian;
fx = fx_current;

  function [x, fx] = hamiltonian(varargin)
  
  x = x_current;
  
  p = normrnd(0,1,size(x_current));
  p_current = p;
  
  lp_current = logpdf_current(varargin{:});
  dlp_current = dlogpdf_current(varargin{:});
  
% $$$   mycheckgrad(@chkgrad, x, 1e-6)
% $$$   
% $$$     function [y,dy] = chkgrad(x)
% $$$     [fx, dfx] = func_x(x);
% $$$     f = get_logpdf(fx);
% $$$     df = get_dlogpdf(dfx);
% $$$     y = f(varargin{:});
% $$$     dy = df(varargin{:});
% $$$     end

  % Random step size
  e = exprnd(epsilon);
  
  % Make a half step for momentum at the beginning
  p = p + 0.5*e*dlp_current;
  
  for l=1:L
    x = x + e * p;
    
    if any(isnan(x))
      break;
    end
    
    [fx, dfx] = func_x(x);
    dlogpdf = get_dlogpdf(dfx);
    dlp = dlogpdf(varargin{:});
    if l < L
      p = p + e * dlp;
    else
      logpdf = get_logpdf(fx);
      lp = logpdf(varargin{:});
    end
  end
  
  if all(~isnan(x))

    % Make a half step for momentum at the end
    p = p + 0.5*e*dlp;
  
    % Negate momentum to make the proposal symmetric
    p = -p;
  
  
    if log(rand()) < ( (lp - 0.5*(p'*p)) - ...
                       (lp_current - 0.5*(p_current'*p_current)) )
      % Accept
      disp('Accept hamiltonian')
      x_current = x;
      fx_current = fx;
      dfx_current = dfx;
      logpdf_current = logpdf;
      dlogpdf_current = dlogpdf;
    else
      disp('Reject hamiltonian')
    end
  
  else
      disp('Reject hamiltonian')
  end
  
  x = x_current;
  fx = fx_current;
  end
  
  
end

function [fx, dfx] = func_x_default(x)
fx = x;
dfx = x;
end

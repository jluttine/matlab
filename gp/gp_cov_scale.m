% GP_COV_SCALE - Adds a scaling parameter to a covariance function.
%
% COVFUNC_SCALED = GP_COV_SCALE(COVFUNC)
%
% [K, DK] = COVFUNC_SCALED(THETA)
%
% Adds the scaling parameter to the beginning of the parameter vector,
% THETA(1), and passes the remaining elements THETA(2:END) to the original
% covariance function.
%
% [N_THETA, N1, N2] = COVFUNC_SCALED()

% Last modified 2010-01-21
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function func = gp_cov_scale(covfunc, varargin)

options = struct('scale', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

func = @get_covariance;

function varargout = get_covariance(theta)
  out = cell(nargout,1);
  varargout = cell(nargout,1);
  
  if nargin == 0
    [out{:}] = covfunc();
    if isempty(options.scale)
      varargout{1} = 1 + out{1};
    else
      varargout{1} = out{1};
    end
    varargout(2:end) = out(2:end);
    return
  end
  
  if isempty(options.scale)
    % Scale parameter in theta
    scale = theta(1);
    n0 = 2;
  else
    % Fixed scale
    scale = options.scale;
    n0 = 1;
  end
  
  [out{:}] = covfunc(theta(n0:end));
  varargout{1} = scale^2 * out{1};
  
  if nargout >= 2
    varargout{2} = cell(numel(theta),1);
    if isempty(options.scale)
      varargout{2}{1} = 2*scale*out{1};
    end
    for n=n0:numel(theta)
      % Scale the gradients
      varargout{2}{n} = scale^2 * out{2}{n-n0+1};
    end
  end
  
  end
  
end

% GP_COV_RQ - Rational quadratic covariance function.
%
% Usage:
%
%   COVFUNC = GP_COV_RQ(D2)
%
% where D2 is a squared distance matrix. The returned covariance function is
% called as
%
%   K = COVFUNC(THETA)
%   [K, DK] = COVFUNC(THETA)
%
% where THETA is a vector of parameters.
%
% THETA(1) : Length scale
% THETA(2) : Alpha, mixture parameter, degrees of freedom
%
% Parameters can be fixed:
%
% [...] = COVFUNC(..., 'LENGTHSCALE', LS)
% [...] = COVFUNC(..., 'ALPHA', ALPHA)
%
% If no arguments are given, dimensionality information is returned:
%
%   [N_THETA, N1, N2] = COVFUNC()
%
% N_THETA : Number of non-fixed parameters
% N1      : Number of rows in the covariance matrix
% N2      : Number of columns in the covariance matrix
%
% See also GP_COV, GP_COV_SE.

% Last modified 2010-01-25
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_rq(D2, varargin)

options = struct('lengthscale', [], ...
                 'alpha', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

n_theta = isempty(options.lengthscale) + isempty(options.alpha);

covfunc = @get_covariance;

  function varargout = get_covariance(theta)

  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      % number of parameters
      varargout{1} = n_theta;
      if nargout >= 2
        varargout{2} = size(D2,1); % number of rows
        if nargout == 3
          varargout{3} = size(D2,2); % number of columns
        else
          error('Too many outputs');
        end
      end
    end
    return
  end

  if numel(theta) ~= n_theta
    error('Wrong number of parameters (%d), should be %d', numel(theta), ...
          n_theta);
  end
  
  % Parse parameters
  if isempty(options.lengthscale)
    ls = theta(1);
    if isempty(options.alpha)
      alpha = theta(2);
    else
      alpha = options.alpha;
    end
  else
    ls = options.lengthscale;
    if isempty(options.alpha)
      alpha = theta(1);
    else
      alpha = options.alpha;
    end
  end

  f = 1 + D2 / (2*alpha*ls^2);
  K = f.^(-alpha);
  varargout{1} = K;

  if nargout >= 2
    dK = cell(n_theta,1);
    n = 1;
    if isempty(options.lengthscale)
      dK{n} =  (K./f.*D2) .* ((-alpha).*(-2)/(2*alpha)*ls^(-3));
% $$$       dK{n} = (-alpha) .* f.^(-alpha-1) .*D2 .* ((-2)/(2*alpha)*ls^(-3));
      n = n + 1;
    end
    if isempty(options.alpha)
      dK{n} = K .* (-log(f)) + ...
              (K./f.*D2) .* ((-1).*(-alpha)./(2*ls^2).*alpha^(-2));
% $$$       dK{n} = K .* (-log(f)) + ...
% $$$               (-alpha) .* f.^(-alpha-1) .* (-1) .* D2 ./ (2*ls^2) .* alpha^(-2);
    end
    varargout{2} = dK;
  end
  
  end
  
end

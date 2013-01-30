% GP_COV_PERIODIC - Squared exponential periodic covariance function.
%
% Usage:
%
%   COVFUNC = GP_COV_PERIODIC(D)
%
% where D is a distance matrix. The returned covariance function is
% called as
%
%   K = COVFUNC(THETA)
%   [K, DK] = COVFUNC(THETA)
%
% where THETA is a vector of parameters.
%
% THETA(1) : Smoothness (length scale relative to the wave length)
% THETA(2) : Wave length
%
% Parameters can be fixed:
%
% [...] = COVFUNC(..., 'SMOOTHNESS', S)
% [...] = COVFUNC(..., 'WAVELENGTH', WL)
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

function covfunc = gp_cov_periodic(D, varargin)

options = struct('smoothness', [], ...
                 'wavelength', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

n_theta = isempty(options.smoothness) + isempty(options.wavelength);

covfunc = @get_covariance;

  function varargout = get_covariance(theta)

  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      % number of parameters
      varargout{1} = n_theta;
      if nargout >= 2
        varargout{2} = size(D,1); % number of rows
        if nargout == 3
          varargout{3} = size(D,2); % number of columns
        else
          error('Too many outputs');
        end
      end
    end
    return
  end

  if numel(theta) ~= 1
    error('Wrong number of parameters (%d), should be %d', numel(theta), ...
          n_theta);
  end
  
  % Parse parameters
  if isempty(options.smoothness)
    ls = theta(1); % smoothness
    if isempty(options.wavelength)
      wl = theta(2); % wave length
    else
      wl = options.wavelength;
    end
  else
    ls = options.smoothness;
    if isempty(options.wavelength)
      wl = theta(1); % wave length
    else
      wl = options.wavelength;
    end
  end

  % Covariance matrix

  angle = pi*D/wl;
  exponent = -2*sin(angle).^2/ls^2;
  K = exp(exponent);
  varargout{1} = K;

  % Gradient for hyperparameters
  if nargout >= 2
    dK = cell(n_theta,1);
    n = 1;
    if isempty(options.smoothness)
      dK{n} = K .* exponent .* (-2) / ls;
      n = n + 1;
    end
    if isempty(options.wavelength)
      dK{n} = K .* (-2/ls^2*2*sin(angle)) .* cos(angle) .* (-pi*D/wl^2);
    end
    varargout{2} = dK;
  end
  end

end

  % TODO:
  % helper, inverse squared length scale:
% $$$   invl2 = theta(1)^(-2); % exp(-2*logtheta(1));
% $$$   wavel = theta(2); %exp(logtheta(2));
  
% $$$ % Gradients for inputs x2
% $$$ if nargout >= 3
% $$$   if isempty(x2)
% $$$     error('Can''t calculate gradient: x2 not given');
% $$$   end
% $$$   d = rows(x2); % dimensionality of inputs
% $$$   m = cols(x1); % number of other inputs
% $$$   n = cols(x2); % number of inputs
% $$$   dK_x2 = zeros([d,m,n]);
% $$$   dK_x2 = bsxfun(@times, ...
% $$$                  reshape(K,[1,m,n]), ...
% $$$                  reshape((-2*invl2*2*sin(angle)) .* ...
% $$$                  cos(angle) .* pi/wavel .* (-1), [1,m,n]));
% $$$ end


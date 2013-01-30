% GP_COV_SE - Squared exponential covariance function.
%
% COVFUNC = GP_COV_SE(D2)
%
% [K, DK] = COVFUNC(THETA)
%
% [N_THETA, N1, N2] = COVFUNC()

% Last modified 2011-01-21
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_se(D2, varargin)

options = struct('lengthscale', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

n_theta = double(isempty(options.lengthscale));

covfunc = @get_covariance;

  function varargout = get_covariance(theta)
  
  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      varargout{1} = n_theta; % number of parameters
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
    error('Wrong number of parameters (%d), should be %d', numel(theta), n_theta);
  end
  
  % Parse parameters
  if isempty(options.lengthscale)
    ls = theta(1);
  else
    ls = options.lengthscale;
  end

  % Covariance matrix
  K = exp(-0.5/ls^2*D2);
  varargout{1} = K;

  % Gradient for hyperparameters
  if nargout >= 2
    dK = cell(n_theta,1);
    if n_theta >= 1
      dK{1} = K .* (ls^(-3)*D2);
    end
    varargout{2} = dK;
  end
  
  end

end

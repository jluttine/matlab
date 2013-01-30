% GP_COV_DAMPEDCOS - Damped cosine covariance function.
%
% c(r) = exp(-r/lengthscale) * cos( 2*pi*r/wavelength )
%
% COVFUNC = GP_COV_DAMPEDCOS(D2)
%
% [K, DK] = COVFUNC(THETA)
%
% [N_THETA, N1, N2] = COVFUNC()

function covfunc = gp_cov_dampedcos(D, varargin)

options = struct('wavelength', [],...
                 'lengthscale', [] );
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

n_theta = isempty(options.wavelength) + isempty(options.lengthscale);

covfunc = @get_covariance;

  function varargout = get_covariance(theta)
  
  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      varargout{1} = n_theta; % number of parameters
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

  if numel(theta) ~= n_theta
    error('Wrong number of parameters (%d), should be %d', numel(theta), n_theta);
  end
  
  % Parse parameters
  if isempty(options.wavelength)
      wl = theta(1);
      if isempty(options.lengthscale)
          ls = theta(2);
      else
          ls = options.lengthscale;
      end
  else
      wl = options.wavelength;
      if isempty(options.lengthscale)
          ls = theta(1);
      else
          ls = options.lengthscale;
      end
  end

  % Covariance matrix
  K = exp(-D/ls) .* cos( 2*pi*D/wl );
  varargout{1} = K;

  % Gradient for hyperparameters
  if nargout >= 2
      dK = cell(n_theta,1);
      n = 1;
      if isempty(options.wavelength)
          dK{n} = exp(-D/ls) .* sin( 2*pi*D/wl ) .* ( 2*pi*D*wl^(-2) );
          n = n + 1;
      end
      if isempty(options.lengthscale)
          dK{n} = exp(-D/ls) .* ( D*ls^(-2) ) .* cos( 2*pi*D/wl );
      end
      varargout{2} = dK;
  end
  
  end % function get_covariance

end

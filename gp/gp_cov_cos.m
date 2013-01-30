% GP_COV_COS - Cosine covariance function.
%
% Is this valid covariance function?
%
% COVFUNC = GP_COV_COS(D)
%
% Length scale can be fixed:
% COVFUNC = GP_COV_COS(D, 'wavelength', 0.7)

% Last modified 2012-01-10
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function covfunc = gp_cov_cos(D, varargin)

options = struct('wavelength', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

n_theta = double(isempty(options.wavelength));

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
  else
    wl = options.wavelength;
  end

  % Covariance matrix
  K = cos( 2*pi*D/wl );
  varargout{1} = K;

  % Gradient for hyperparameters
  if nargout >= 2
    dK = cell(n_theta,1);
    if n_theta >= 1
      dK{1} = sin( 2*pi*D/wl ) .* ( 2*pi*D*wl^(-2) );
    end
    varargout{2} = dK;
  end
  
  end

end

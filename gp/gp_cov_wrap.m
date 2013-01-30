% GP_COV_WRAP - A covariance function which wraps a given covariance
%               matrix as a covariance function.
%
% The function makes it easy to create a dummy covariance function wrapper
% for a covariance matrix. The covariance function simply returns the given
% covariance matrix.
%
% COVFUNC = GP_COV_WRAP(C)
%
% [K, DK] = COVFUNC(THETA) 
%
% The covariance function COVFUNC has no parameters (THETA is empty) and it
% returns K=C. The derivative DK is an empty cell.
%
% [N_THETA, N1, N2] = COVFUNC()

% Last modified 2010-01-21
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_wrap(C)

covfunc = @get_covariance;

  function varargout = get_covariance(theta)
  
  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      varargout{1} = 0; % number of parameters
      if nargout >= 2
        varargout{2} = size(C,1); % number of rows
        if nargout == 3
          varargout{3} = size(C,2); % number of columns
        else
          error('Too many outputs');
        end
      end
    end
    return
  end

  % Covariance matrix
  varargout{1} = C;

  % Gradient for hyperparameters
  if nargout >= 2
    varargout{2} = cell(0,1);
  end
  
  end
  
end
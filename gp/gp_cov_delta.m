
function covfunc = gp_cov_delta(N)

covfunc = @get_covariance;

  function varargout = get_covariance(theta)
  
  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      varargout{1} = 0; % number of parameters
      if nargout >= 2
        varargout{2} = N; % number of rows
        if nargout == 3
          varargout{3} = N; % number of columns
        else
          error('Too many outputs');
        end
      end
    end
    return
  end

  K = speye(N);
  varargout{1} = K;
  if nargout >= 2
    dK = cell(0,1);
    varargout{2} = dK;
  end
  end

end
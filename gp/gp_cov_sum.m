% GP_COV_SUM - Covariance function as a sum of two covariance functions.
%
% Usage:
%
%   COVFUNC = GP_COV_SUM(COVFUNC1, COVFUNC2)
%
% The returned covariance function is called as
%
%   K = COVFUNC(THETA)
%   [K, DK] = COVFUNC(THETA)
%
% where THETA is a vector of parameters. The parameters of COVFUNC1 and
% COVFUNC2 are concatenated into a single parameter vector.
%
% If no arguments are given, dimensionality information is returned:
%
%   [N_THETA, N1, N2] = COVFUNC()
%
% N_THETA : Number of non-fixed parameters
% N1      : Number of rows in the covariance matrix
% N2      : Number of columns in the covariance matrix
%
% See also GP_COV, GP_COV_PRODUCT.

% Last modified 2010-01-25
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function func = gp_cov_sum(varargin)

func = @get_covariance;

if nargin < 2
  error('Must give at least two covariance functions')
end

covfuncs = varargin ;

% Number of covariance functions
n_funcs = nargin ;

% Number of parameters for each covariance function
n_theta = zeros(n_funcs,1) ;
% Dimensionalities of each covariance function
M = zeros(n_funcs,1) ;
N = zeros(n_funcs,1) ;
% Extract the values from the covariance functions
for i = 1:n_funcs
  [n_theta(i), M(i), N(i)] = varargin{i}();
end

% Indices of the hyperparameters for each covariance function
ind_theta = cell(n_funcs,1);
for i = 1:n_funcs
  ind_theta{i} = (1+sum(n_theta(1:(i-1)))):(sum(n_theta(1:i))) ;
end

%[n_theta1, M1, N1] = covfunc1();
%[n_theta2, M2, N2] = covfunc2();

if ~all(M == M(1)) || ~all(N == N(1))
  error('Can''t sum covariance matrices with different dimensionalities');
end
M = M(1) ;
N = N(1) ;

  function varargout = get_covariance(theta)
  
  if nargout == 0
    nout = 1 ;
  else
    nout = nargout ;
  end
  
  varargout = cell(nout,1);

  out = cell(nout,1);
  varargout = cell(nout,1);
  
  % Return only dimension information if requested
  if nargin == 0
    if nout >= 1
      varargout{1} = sum(n_theta); % number of parameters
      if nout >= 2
        varargout{2} = M; % dimensionalities
        if nout >= 3
          varargout{3} = N; % dimensionalities
        end
      end
    end
    return
  end

  varargout{1} = 0 ;
  varargout{2} = cell(sum(n_theta),1) ;
  
  for i = 1:n_funcs
    % Compute the covariance matrix
    [out{:}] = covfuncs{i}(theta(ind_theta{i}));
    varargout{1} = varargout{1} + out{1};
  
    % Compute the derivative
    if nout >= 2
      varargout{2}(ind_theta{i}) = out{2}(:);
    end
  end
  
  end

end
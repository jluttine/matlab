% GP_COV_KRON - Covariance function as a product of two covariance
%                  functions.
%
% Usage:
%
%   COVFUNC = GP_COV_PRODUCT(COVFUNC1, COVFUNC2, ...)
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
% See also GP_COV, GP_COV_SUM.

% Last modified 2012-03-14
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function func = gp_cov_kron(varargin)

if nargin < 2
  error('Must give at least two covariance functions')
end

func = @get_covariance;

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

%[n_theta2, M2, N2] = covfunc2();

%if M1 ~= M2 || N1 ~= N2
%  error('Can''t multiply covariance matrices with different dimensionalities');
%end

  function varargout = get_covariance(theta)
  
  if nargout == 0 
    nout = 1;
  else
    nout = nargout;
  end
  
  varargout = cell(nout,1);

  out = cell(nout,1);
  %out2 = cell(nout,1);
  varargout = cell(nout,1);
  
  % Return only dimension information if requested
  if nargin == 0
    if nout >= 1
      varargout{1} = sum(n_theta); % number of parameters
      if nout >= 2
        varargout{2} = prod(M) ; % dimensionalities
        if nout >= 3
          varargout{3} = prod(N); % dimensionalities
        end
      end
    end
    return
  end

  if numel(theta) ~= sum(n_theta)
    error('Wrong number of parameters (%d), should be %d', numel(theta), ...
          sum(n_theta));
  end
  
  % Initialize covariance matrix and gradients
  varargout{1} = 1 ;
  if nout >= 2
    varargout{2} = cell(sum(n_theta),1) ;
    for i = 1:sum(n_theta)
      varargout{2}{i} = 1 ;
    end
  end
  

  % Compute the covariance matrix
  for i = 1:n_funcs
    [out{:}] = covfuncs{i}(theta(ind_theta{i}));
    varargout{1} = kron(varargout{1}, out{1});
  
    % Compute the derivative
    if nout >= 2 && ~isempty(ind_theta{i})
      % No derivative for these hyperparameters
      for n = 1:(ind_theta{i}(1)-1)
        varargout{2}{n} = kron(varargout{2}{n}, out{1}) ;
      end
      % Derivative for these hyperparameters
      for n = 1:n_theta(i)
        ind = ind_theta{i}(1) + n - 1 ;
        varargout{2}{ind} = kron(varargout{2}{ind}, out{2}{n}) ;
      end
      % No derivative for these hyperparameters
      for n = (ind_theta{i}(end)+1):sum(n_theta)
        varargout{2}{n} = kron(varargout{2}{n}, out{1}) ;
      end
    end
  end
  
  end


end

%function X = kronecker(varargin)

%X = varargin{1} ;
%for i = 2:nargin
%  X = kron(X, varargin{i}) ;
%end

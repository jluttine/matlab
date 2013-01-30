% GP_COV_PRODUCT - Covariance function as a product of two covariance
%                  functions.
%
% Usage:
%
%   COVFUNC = GP_COV_PRODUCT(COVFUNC1, COVFUNC2)
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

% Last modified 2010-01-25
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_product(covfunc1, covfunc2)

covfunc = @get_covariance;

[n_theta1, M1, N1] = covfunc1();
[n_theta2, M2, N2] = covfunc2();

if M1 ~= M2 || N1 ~= N2
  error('Can''t multiply covariance matrices with different dimensionalities');
end

  function varargout = get_covariance(theta)
  
  varargout = cell(nargout,1);

  out1 = cell(nargout,1);
  out2 = cell(nargout,1);
  varargout = cell(nargout,1);
  
  % Return only dimension information if requested
  if nargin == 0
    if nargout >= 1
      varargout{1} = n_theta1 + n_theta2; % number of parameters
      if nargout >= 2
        varargout{2} = M1; % dimensionalities
        if nargout >= 3
          varargout{3} = N1; % dimensionalities
        end
      end
    end
    return
  end

  % Compute the covariance matrix
  [out1{:}] = covfunc1(theta(1:n_theta1));
  [out2{:}] = covfunc2(theta((n_theta1+1):end));
  varargout{1} = out1{1} .* out2{1};
  
  % Compute the derivative
  if nargout >= 2
    dK = cell(numel(theta),1);
    for n=1:n_theta1
      dK{n} = out1{2}{n} .* out2{1};
    end
    for n=1:n_theta2
      dK{n+n_theta1} = out1{1} .* out2{2}{n};
    end
    varargout{2} = dK;
  end
  end

end
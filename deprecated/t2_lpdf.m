% Student-t log-probability density function (lpdf) with degrees of
% freedom NU.
%
% Y = T2_LPDF(X2, NU, D)
% [Y, DY] = T2_LPDF(X2, NU, D)
%
% DY is the derivative with respect to the degrees of freedom NU.
%
% X2 is the squared values of a Student-t random variable with zero mean
% and unit scale.
%
% or in the multivariate case:
%
% ... = T2_LPDF(..., D)
%
% where D is the dimensionality, and X2 is the squared norm of a
% Student-t random variable with zero mean and unit scale.

% Last modified 2010-06-03
% Copyright (c) Jaakko Luttinen

% X2 = ((T-MU)/SIGMA)^2 is the normalized squared error.
% X2 = (T-MU)' * inv(COV) * (T-MU)

function [y, dy] = t2_lpdf(x2, nu, D)
warning('This function is deprecated')


warning('Deprecated. Use T_LOGPDF instead.')

if nargin < 3
  D = 1;
end

b = (nu+D)/2;
z = log(1 + 1./nu.*x2);

y = gammaln(b) - gammaln(nu/2) - D./2.*log(nu*pi) - b .* z;

if nargout >= 2
  dy = 0.5*psi(b) ...
       - 0.5*psi(nu/2) ...
       - 0.5*D./nu ...
       - 0.5*z ...
       + b.*x2./(nu.*(nu+x2));
% $$$   dy = 0.5 * ( psi(b) - psi(nu/2) - D./nu - z + (nu+1).*x2./ (nu.*(nu+x2)) );
end


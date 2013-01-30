function [y, dy] = t2_lpdf(x2, nu)
% Student-t log-probability density function (lpdf) with degrees of
% freedom NU.
%
% Y = T2_LPDF(X2, NU)
% [Y, DY] = T2_LPDF(X2, NU)
%
% DY is the derivative with respect to the degrees of freedom NU.
%
% X2 = ((T-MU)/SIGMA)^2 is the normalized squared error.

error('Not in use?')

b = (nu+1)/2;
z = log(1 + 1./nu.*x2);

y = gammaln(b) - gammaln(nu/2) - 0.5*log(nu*pi) - (nu+1)/2 .* z;

disp('debuggin too')
if nargout >= 2
  disp('debuggin')
  dy = 0.5 * ( psi(b) - psi(nu/2) - 1./nu - z + (nu+1).*x2./ (nu.*(nu+x2)) );
end


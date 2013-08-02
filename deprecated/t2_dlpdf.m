function dy = t2_dlpdf(x2, nu)
% Derivative of the Student-t log-probability density function (lpdf) with
% degrees of freedom NU.
%
% DY = T2_LPDF(X2, NU)
%
% DY is the derivative with respect to the degrees of freedom NU.
%
% X2 = ((T-MU)/SIGMA)^2 is the normalized squared error.

b = (nu+1)/2;
z = log(1 + 1./nu.*x2);

dy = 0.5 * ( psi(b) - psi(nu/2) - 1./nu - z + (nu+1).*x2./ (nu.*(nu+x2)) );



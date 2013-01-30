% T_LOGPDF - The log pdf of a multivariate Student-t distribution.
%
% Computes the logarithm of the (multivariate) Student-t probability density
% function STUDENT-T(Y|MU,COV,NU), where Y is the random variable, MU is the
% location, COV is the scale matrix and NU is the degrees of freedom.
%
% The function is called as:
%
% L = T_LOGPDF(X2, LOGDET_COV, NU, D)
% [L, DL] = T_LOGPDF(X2, LOGDET_COV, NU, D)
%
% where
%
% X2         : (Y-MU)'*INV(COV)*(Y-MU)
% LOGDET_COV : LOG(DET(COV))
% NU         : the degrees of freedom
% D          : the dimensionality of the distribution (default is 1)
%
% L  : LOG(STUDENT-T(Y|MU,COV,NU))
% DL : the derivative of L with respect to NU (optional)

% Last modified 2010-06-08
% Copyright (c) Jaakko Luttinen

function [l, dl] = t_logpdf(x2, logdet_Cov, nu, D)

if nargin < 3
  D = 1;
end

lognu = log(nu);
logx2 = log(x2);

b = (nu+D)/2;
z = log(1 + exp(-lognu+logx2));
% z = log(1 + 1./nu.*x2);

l = gammaln(b) - gammaln(nu/2) - 0.5*logdet_Cov - D./2.*(lognu+log(pi)) - b .* z;

if nargout >= 2
  dl = 0.5*psi(b) ...
       - 0.5*psi(nu/2) ...
       - 0.5*D./nu ...
       - 0.5*z ...
       + exp( log(b) + logx2 - (lognu + log(nu+x2)) );
%       + b.*x2./(nu.*(nu+x2));
end


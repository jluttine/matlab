function c = mnorm_vbcost(x, L, mu, L_mu, L_K, logdet_K)
warning('This function is deprecated')

% C = MNORM_VBCOST(X, L, MU, L_MU, L_K, LOGDET_K)
%
% Evaluates C = -KL(q||p) = <log p(X)> - <log q(X)> where the expectation
% <.> is taken over q. This is used for calculating the terms of the lower
% bound of the log marginal likelihood in variational Bayesian inference.
%
% p(X) = N(M, K)
% q(X) = N(X, L*L')
% 
% Rest of the parameters are defined as:
% q(M) = N(MU, L_MU * L_MU')
% <inv K> = inv(L_K * L_K')
% <log det K> = LOGDET_K
%
% If L==[], then the term <log q(X)> is not calculated.
% This is useful, for instance, when X is observations.

% (c) 2010 Jaakko Luttinen


lq = 0;
d = length(x);

%% Cost from approximate posterior

if ~isempty(L)
  lq = mnorm_entropy(L);
else
  L = 0;
end

%% Cost from prior

% Use Cholesky of the prior covariance
z = linsolve_tril(L_K, x-mu);
V = linsolve_tril(L_K, L+L_mu)
z2 = z'*z + traceprod(V, V);

lp = - 0.5*logdet_K - 0.5*z2 - 0.5*d*log(2*pi);

%% Total cost
c = lp + lq;
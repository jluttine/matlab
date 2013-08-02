function c = mvn_vbcost(x,L,mu,L_mu,L_K,logdet_K)
warning('This function is deprecated')

% c = mvnvbcost(x,Covx,mu,Covmu,pCholCov,pLogDetCov)
%
% Calculates c = -KL(q||p) = <log p(X)> - <log q(X)>
% where the expectation <.> is taken over q(X). This is used for
% calculating the lower bound of the marginal loglikelihood in VB
% models.
%
% p(X) = N(M,K)
% q(X) = N(x,Covx)
% 
% Rest of the parameters are defined as:
% q(M) = N(mu,Covmu)
% <inv K> = inv(L_K * L_K')
% <log det K> = logdet_K
%
% If Covx==[], then the term <log q(X)> is not calculated.
% This is useful when X is, e.g., observations.


c = 0;
d = length(x);

% Cost from q-posterior
if ~isempty(L)
  c = mnorm_entropy(L);
else
  L = 0;
end

% Use Cholesky of the prior covariance
z = solve_tril(L_K, x-mu);

% Cost from prior

% Below: (x-mu)' * Cov^(-1) * (x-mu) + trace(Cov^(-1)*(Covx+Covmu))
V = linsolve_chol(L_K, L*L'+L_mu*L_mu')
err2 = z'*z + trace(V, V); % TODO: Optimize!
% $$$ err2 = z'*z + trace(solve_triu(pCholCov',solve_tril(pCholCov,Covx+Covmu)));
c = c - 0.5*logdet_K - 0.5*err2 - 0.5*d*log(2*pi);

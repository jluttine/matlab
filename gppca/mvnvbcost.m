function c = mvnvbcost(x,Covx,mu,Covmu,pCholCov,pLogDetCov)
% c = mvnvbcost(x,Covx,mu,Covmu,pCholCov,pLogDetCov)
%
% Calculates c = -KL(q||p) = <log p(X)> - <log q(X)>
% where the expectation <.> is taken over q(X). This is used for
% calculating the lower bound of the marginal loglikelihood in VB
% models.
%
% p(X) = N(M,T)
% q(X) = N(x,Covx)
% 
% Rest of the parameters are defined as:
% q(M) = N(mu,Covmu)
% <inv T> = inv(pCholCov * pCholCov')
% <log det T> = pLogDetCov
%
% If Covx==[], then the term <log q(X)> is not calculated.
% This is useful when X is, e.g., observations.


c = 0;
d = length(x);

% Cost from q-posterior
if ~isempty(Covx)
  % TODO:
  % if ~isempty(pCholCov), entropy = mvnentropy(pCholCov,true); else.... end
%  entropy = mvnentropy(pCholCov,true);
  entropy = mvnentropy(Covx);
  c = entropy;
else
  Covx = 0;
end

% Use Cholesky of the prior covariance
% $$$ opts.UT = true;
% $$$ opts.TRANSA = true;
z = solve_tril(pCholCov, x-mu);

% Cost from prior

% Below: (x-mu)' * Cov^(-1) * (x-mu) + trace(Cov^(-1)*(Covx+Covmu))
err2 = z'*z + trace(solve_triu(pCholCov',solve_tril(pCholCov,Covx+Covmu)));
%err2 = (x-mu)'*pInvCov*(x-mu) + (Covx(:)+Covmu(:))'*pInvCov(:);
%((x-mu).^2 + diag(Covx) + diag(Covmu));
c = c - 0.5*pLogDetCov - 0.5*err2 - 0.5*d*log(2*pi);

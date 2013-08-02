
function c = gam_vbcost(x,logx,apost,bpost,aprior,bprior)
% function c = gam_vbcost(x,logx,apost,bpost,aprior,bprior)
%
% Calculates c = -KL(q||p) = <log p(X)> - <log q(X)>
% where the expectation <.> is taken over q(X). This is used for
% calculating the lower bound of the marginal loglikelihood in VB
% models.
%
% Let X be a Gamma distributed variable:
% p(X) = G(aprior,bprior)
% q(X) = G(apost,bpost)
% <X> = x
% <log X> = logx

% Cost from prior
c = aprior.*log(bprior) - gammaln(aprior) + (aprior-1).*logx - bprior.*x;
% Cost from q-posterior
c = c + gam_entropy(apost,bpost);

function d = gamkl(x,logx,apost,bpost,aprior,bprior)
warning('This function is deprecated')

% KL divergence for Gamma distributions.
% A lower bound for a Gamma distributed X:
% p(X) = G(aprior,bprior)
% q(X) = G(apost,bpost)
% <X> = x
% <log X> = logx
%
% KL(q||p) = <log q(X)> - <log p(X)> = -H(q) - <log p(x)>
% where expectation is taken over q(X) and H() is entropy

% Cost from prior
d = -aprior*log(bprior) + gammaln(aprior) - (aprior-1)*logx + bprior*x;
% Cost from q-posterior
d = d - gament(apost,bpost);

function p = mstud_lpdf(x, mu, Cov, nu)
% p = mstud_lpdf(x, mu, Cov, nu)

D = length(x);

x = x(:);
mu = mu(:);

invCov = inv(Cov);

p = gammaln((D+nu)/2);
p = p - gammaln(nu/2);
p = p - D/2 * log(nu * pi);
p = p - 0.5 * log(det(Cov));
p = p + -(D+nu)/2 * log(1 + 1/nu*(x-mu)'*invCov*(x-mu));
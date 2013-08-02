
% function y = mnormrnd(mu, Cov, n)
function y = mnormrnd(mu, Cov, n)

m = max(length(mu), length(Cov));

opts.issym = true;
opts.isreal = true;
[V,D] = svd(Cov);
D(D<0) = 0;
D = sqrt(D);
A = V * D;

y = repmat(mu,1,n) + A*randn(m,n);
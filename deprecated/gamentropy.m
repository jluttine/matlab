function H = gamentropy(a,b)
warning('This function is deprecated')

% H = entropy_gamma(a,b)
% a is shape, b is inverse scale, that is, the mean is a/b.

H = a - log(b) + gammaln(a) - (a-1).*psi(a);


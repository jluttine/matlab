function H = gament(a,b)
% a is shape, b is inverse scale
% Entropy of a Gamma distribution.

H = a - log(b) + gammaln(a) - (a-1).*psi(a);


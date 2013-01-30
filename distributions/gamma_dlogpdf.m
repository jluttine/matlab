% GAMMA_DLOGPDF - The derivative of the log pdf of a Gamma distribution.
%
% Computes the derivative of the logarithm of the Gamma probability density
% function G(ALPHA|A,B), where ALPHA is the random variable, A is the shape
% and B is the inverse scale (note the difference to MATLAB's built-in
% functions).
%
% The function is called as
%
% D = GAMMA_DLOGPDF(ALPHA, A, B)
% D = GAMMA_DLOGPDF(DALPHA, A, B, DLOG_ALPHA)

% Last modified 2011-02-09
% Copyright (c) Jaakko Luttinen

function dlp = gamma_dlogpdf(alpha, a, b, dlog_alpha)
%function dlp = gamma_dlogpdf(alpha, a, b)

if nargin < 4
  inv_alpha = 1./alpha;
  dlp = (a-1).*inv_alpha - b;
else
  dlp = (a-1).*dlog_alpha - b.*alpha;
end

% GAMMA_LOGPDF - The log pdf of a Gamma distribution.
%
% Computes the logarithm of the Gamma probability density function
% G(ALPHA|A,B), where ALPHA is the random variable, A is the shape and B is
% the inverse scale (note the difference to MATLAB's built-in functions).
%
% The function is called as
%
% L = GAMMA_LOGPDF(ALPHA, A, B)
% L = GAMMA_LOGPDF(ALPHA, A, B, LOG_ALPHA)
%
% where LOG_ALPHA is LOG(ALPHA). Letting the user to evaluate this logarithm
% allows greater flexibility (e.g., for variational Bayesian cost function
% one can give the expectation of the logarithm). However, A and B are
% assumed to be fixed parameters in order to keep the function simple.

% Last modified 2010-06-07
% Copyright (c) Jaakko Luttinen

function [lp, dlp] = gamma_logpdf(alpha, a, b, log_alpha)

if nargin < 4
  log_alpha = log(alpha);
% $$$   if nargin == 4
% $$$     error(['Check the order of your input variables - the function has been ' ...
% $$$            'changed!'])
% $$$   end
% $$$ else
% $$$   error('Improper number of inputs')
end

lp = (a-1).*log_alpha - b.*alpha - gammaln(a) + a.*log(b);

if nargout >= 2
  inv_alpha = 1./alpha;
  dlp = (a-1).*inv_alpha - b;
end

%lp(alpha<=0) = -inf;
%dlp(alpha<=0) = nan;
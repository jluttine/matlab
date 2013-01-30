% NORMAL_LOGPDF
%
% Y = NORMAL_LOGPDF(X)
% Y = NORMAL_LOGPDF(X, MU, SIGMA)
%
% MU is the location parameter in (-INF, INF), default 0
% SIGMA is the scale parameter in (0, INF)

% Last modified 2010-11-12
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function y = lognormal_logpdf(x, varargin)

if nargin <= 3
  if nargin == 1
    mu = 0;
    sigma = 1;
  else
    mu = varargin{1};
    sigma = varargin{2};
  end
  y = -0.5*log(2*pi) ...
      -log(sigma) ...
      -0.5*((x-mu)./sigma).^2;
end

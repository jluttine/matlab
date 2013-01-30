% LOGNORMAL_RAND
%
%   X = LOGNORMAL_RAND(MU, SIGMA)
%
% MU is the location parameter in (-INF, INF)
% SIGMA is the scale parameter in (0, INF)
%
% Function calls
%   X = LOGNORMAL_RAND(MU,SIGMA,M,N,...)
%   X = LOGNORMAL_RAND(MU,SIGMA,[M,N,...])  
% return an M-by-N-by-... array.
%
% See also LOGNORMAL_LOGPDF, NORMRND.

% Last modified 2010-11-12
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = lognormal_rand(varargin)

x = exp(normrnd(varargin{:}));
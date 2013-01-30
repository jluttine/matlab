% RMSE - Evaluates root mean squared error.
%
% E = RMSE(Y)
% 
% E = RMSE(Y1,Y2)

% Last modified 2010-10-21
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function e = rmse(Y, Yh)

if nargin < 2
  e = sqrt(Y(:)'*Y(:)/numel(Y));
else
  z = Y(:)-Yh(:);
  e = sqrt(z'*z/numel(z));
end


% K = gpcovConstScale(s, covfunc, ...)
% [K, dK_logtheta] = gpcovConstScale(s, covfunc, ...)
% [K, dK_logtheta, dK_x2] = gpcovConstScale(s, covfunc, ...)
%
% For instance,
% covfunc = {@gpcovConstScale, 0.5, {@gpcovSEiso}}
% and then
% K = feval(covfunc{:}, logtheta, x1, x2)

function [K, dK_logtheta, dK_x2] = gpcovNoise(logtheta,x1,x2)

warning('Not well tested / debugged yet');
  
n1 = cols(x1);
if nargin < 3
  % Asking only vector of variances
  n2 = 1;
  s2 = exp(logtheta(1));
  K = s2 * ones(n1,1);
else
  % Asking covariance matrix
  n2 = cols(x2);
  if all(size(x1) == size(x2)) && all(x1 == x2)
    % Noise variance only..
    s2 = exp(logtheta(1));
  else
    % .. no correlations between different inputs
    s2 = 0;
  end
  K = s2 * eye(n1,n2);
end


% Gradient for hyperparameters
if nargout >= 2
  dK_logtheta = zeros([n1, n2, 1]);
  dK_logtheta(:,:,1) = K;
end

% Gradients for inputs x2
if nargout >= 3
  if isempty(x2)
    error('Can''t calculate gradient: x2 not given');
  end
  d = rows(x2); % dimensionality of inputs
  dK_x2 = zeros([d,n1,n2]);
end


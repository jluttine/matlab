
% [K, dK_logtheta, dK_x2] = gpcov(x1, x2, logtheta)
%
function [K, dK_logtheta, dK_x2] = gpcovDot(logtheta, x1, x2)

warning('Under development, not ready')

% Covariance matrix
K = x1' * x2;

% Gradient for hyperparameters
if nargout >= 2
  dK_logtheta = nan;
end

% Gradients for inputs x2
if nargout >= 3
  if isempty(x2)
    error('Can''t calculate gradient: x2 not given');
  end
  d = rows(x2); % dimensionality of inputs
  m = cols(x1); % number of other inputs
  n = cols(x2); % number of inputs
  dK_x2 = zeros([d,m,n]);
  for j=1:n
    dK_x2(:,:,j) = x1;
  end
end


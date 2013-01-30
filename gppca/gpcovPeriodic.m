function [K, dK_logtheta, dK_x2] = gpcovPeriodic(logtheta, x1, x2, distfun)

if nargin == 0
  % Return the number of hyperparameters
  K = 2;
  return
end

if nargin < 3
  x2 = [];
end
if nargin < 4
%  distfun = @sqdistEuclidean;
end

% 2-norm distances
if ~isempty(x2)
  D = bsxfun(@minus, x1(:), x2(:)');
  %dD_dx2 = reshape(-D, [1,length(x1),length(x2)]);
%  [D, dD_dx2] = feval(distfun, x1, x2, false);
else
  % Distances for variances
  n = cols(x1);
  D = zeros(n,1);
end

% TODO:
% helper, inverse squared length scale:
invl2 = exp(-2*logtheta(1));
wavel = exp(logtheta(2));

% Covariance matrix

angle = pi*D/wavel;
exponent = -2*sin(angle).^2*invl2;
K = exp(exponent);

% Gradient for hyperparameters
if nargout >= 2
  dK_logtheta = zeros([size(D), 2]);
  dK_logtheta(:,:,1) = K .* exponent .* (-2);
  dK_logtheta(:,:,2) = K.*(-2*invl2*2*sin(angle)).*cos(angle).*(-pi*D/wavel);
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
  dK_x2 = bsxfun(@times, ...
                 reshape(K,[1,m,n]), ...
                 reshape((-2*invl2*2*sin(angle)) .* ...
                 cos(angle) .* pi/wavel .* (-1), [1,m,n]));
% $$$   for j=1:n
% $$$     dK_x2(:,:,j) = bsxfun(@times, K(:,j)', (-0.5*invl2) * dD2_dx2(:,:,j));
% $$$     %2*bsxfun(@minus, x2(:,j),x1)); 
% $$$   end
end


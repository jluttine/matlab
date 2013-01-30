
% [K, dK_logtheta, dK_x2] = gpcov(logtheta, x1, x2, distfun)
%
function [K, dK_logtheta, dK_x2] = gpcov(logtheta, x1, x2, distfun)

if nargin == 0
  % Return the number of hyperparameters
  K = 1;
  return
end

if nargin < 3
  x2 = [];
end
if nargin < 4
  distfun = @sqdistEuclidean;
end

% 2-norm distances
if ~isempty(x2)
  [D2, dD2_dx2] = feval(distfun, x1, x2, true);
  % sq_dist(x1,x2);
% $$$   % Distances for covariances
% $$$   n1 = cols(x1);
% $$$   n2 = cols(x2);
% $$$   D2 = zeros(n1,n2);
% $$$   % TODO: Could check whether n2 or n1 is smaller and loop over that
% $$$   for j=1:n2
% $$$     Dif = bsxfun(@minus,x1, x2(:,j));
% $$$     D2(:,j) = dot(Dif,Dif,1);
% $$$ %    D2(:,j) = rowsum(bsxfun(@minus,x1, x2(:,j)).^2);
% $$$   end
else
  % Distances for variances
  n = cols(x1);
  D2 = zeros(n,1);
end

% TODO:
% helper, inverse squared length scale:
invl2 = exp(-2*logtheta(1));

% Covariance matrix
K = exp(-0.5*invl2*D2);
% $$$ figure
% $$$ imagesc(K)

% Gradient for hyperparameters
if nargout >= 2
  dK_logtheta = zeros([size(D2), 1]);
  dK_logtheta(:,:,1) = K .* (-0.5*D2*invl2) .* (-2);
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
  dK_x2 = bsxfun(@times, reshape(K,[1,m,n]), (-0.5*invl2) * dD2_dx2);
% $$$   for j=1:n
% $$$     dK_x2(:,:,j) = bsxfun(@times, K(:,j)', (-0.5*invl2) * dD2_dx2(:,:,j));
% $$$     %2*bsxfun(@minus, x2(:,j),x1)); 
% $$$   end
end


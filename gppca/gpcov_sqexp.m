
% [K, dK_logtheta, dK_x2] = gpcov_sqexp(x1, x2, logtheta)
%
function [K, dK_logtheta, dK_x2] = gpcov_sqexp(x1, x2, logtheta)

% 2-norm distances
if ~isempty(x2)
  % Distances for covariances
  n1 = cols(x1);
  n2 = cols(x2);
  D = zeros(n1,n2);
  for i=1:n1
    for j=1:n2
      D(i,j) = norm(x1(:,i) - x2(:,j));
    end
  end
else
  % Distances for variances
  n = rows(x1);
  D = zeros(n,1);
end
  

% Covariance matrix
D2 = D.^2;
K = exp(logtheta(1)) * exp(-0.5*exp(-2*logtheta(2))*D2);

% Gradient for hyperparameters
if nargout >= 2
  dK_logtheta = zeros([size(D), 2]);
  dK_logtheta(:,:,1) = K;
  dK_logtheta(:,:,2) = K .* (-0.5*D2*exp(-2*logtheta(2))) .* (-2);
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
  for i=1:m
    for j=1:n
      dK_x2(:,i,j) = K(i,j) * (-0.5*exp(-2*logtheta(2))) * 2*(x2(:,j)-x1(:,i));
    end
  end
end



function [K, dK] = gpK_sqexp(D, p1, p2)

% [K, dK] = gpK_sqexp(D, p1, p2)
%
% Squared exponential covariance function for GP.
%
% p1 is the scale
% p2 is the length
%
K = p1^2*exp(-0.5*D.^2/(p2^2));

if nargout >= 2
  dK = zeros([size(D), 2]);
  dK(:,:,1) = K .* 2 / p1;
  dK(:,:,2) = K .* (-0.5*D.^2) .* (-2*p2^(-3));
end

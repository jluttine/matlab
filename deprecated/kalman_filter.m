% KALMAN_FILTER - Under construction..

function [X, CovX] = kalman_filter(measurement, ...
                                   transition_matrix, ...
                                   process_noise, ...
                                   measurement_model_matrix, ...
                                   measurement_noise, ...
                                   x, ...
                                   C)

error('Not ready yet..')
return

N = 0;
D = 1000;

X = zeros(D,N);
if nargout >= 2
  CovX = zeros(D,D,10);
end

for n=1:N
  % Get the parameters
  A = transition_matrix(n);
  Q = process_noise(n);
  y = measurement(n);
  H = measurement_model_matrix(n);
  R = measurement_noise(n);
  
  % Prediction step
  x = A*x;
  C = A*C*A' + Q;
  
  % Update step
  if length(y) < D
    v = y - H*x;
    S = H*C*H' + R;
    K = C*H'/S;
    x = m + K*v;
    C = C - K*S*K';
  else
    K = C*H';
    L = chol(C + 
  end
  
  % Store the distribution
  X(:,n) = x;
  if nargout >= 2
    CovX(:,:,n) = C;
  end
end
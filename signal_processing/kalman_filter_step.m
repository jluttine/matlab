
% p(y_n) = N(y_n | C*x_n, R)
% p(x_n | x_(n-1)) = N(x_n | A*x_(n-1), Q)
%
% [x, Covx, logp] = kalman_filter_step(x, Covx, y, A, Q, C, R, 1)
%
% [x, Covx, logp] = kalman_filter_step(x, Covx, y, A, Q, C, inv(chol(R,'lower')), 2)

% Last modified 2011-10-19
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function [x, Covx, logp] = kalman_filter_step(x, Covx, y, A, Q, C, R, method)

% Prediction step
x = A*x;
Covx = A*Covx*A' + Q;

% Update step
switch method
 case 1
  
  % Standard update step (if length(y) <= length(x))
  y = y - C*x;
  S = C*Covx*C' + R;
  L = chol(S, 'lower');
  K = linsolve_tril(L, C*Covx); %(Covx * C') / S;
  v = linsolve_tril(L, y);
  if nargout >= 3
    logp = gaussian_logpdf(v'*v, ...
                           0, ...
                           0, ...
                           logdet_chol(L), ...
                           length(y));
  end
  x = x + K'*v;
  Covx = Covx - K'*K; % K*C*Covx;

 case 2
  % Update step in lower-dimensional space of X
  
  invL_R = R;
  Lx = chol(Covx, 'lower');
  Z = invL_R * C * Lx;
  L = chol(speye(size(Lx)) + Z'*Z, 'lower');
  v_y = y - C*x;
  v_x = C'*(invL_R'*(invL_R*y)) + linsolve_lchol(Lx,x);
  Z = linsolve_tril(L, Lx');
  Covx = Z'*Z;
  x = Covx * v_x;
  
  if nargout >= 3
    v = invL_R' * (invL_R*v_y); % R \ (y - C*x)
    logp = gaussian_logpdf(v'*v_y - v'*C*Covx*C'*v, ...
                           0, ...
                           0, ...
                           logdet_chol(L) - logdet_chol(invL_R), ...
                           length(y));
  end
  
  
% $$$   invR = R;
% $$$   S = Covx * C';
% $$$   K = Covx + S*invR*S';
% $$$   L = chol(K, 'lower');
% $$$   Z = L \ Covx;
% $$$   if nargout >= 3
% $$$     % Compute marginal likelihood term
% $$$   end
% $$$   x = Z' * linsolve_tril(L, S*(invR*y) + x);
% $$$   Covx = Z'*Z;
  
end

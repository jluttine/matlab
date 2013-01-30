% [x, Covx, Covx_x, entropy] = rts_smoother_step(x, Covx, x_s, Covx_s, A, Q)

% Last modified 2011-10-19
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function [x, Covx, Covx_x, entropy] = rts_smoother_step(x, Covx, x_s, Covx_s, A, Q)
% Perform the RTS smoothing step

x_p = A*x;
Covx_p = A*Covx*A' + Q;

S = (Covx*A') / Covx_p;
x = x + S*(x_s-x_p);
if nargout >= 2
  Covx = Covx + S*(Covx_s-Covx_p)*S';
end
if nargout >= 3
  Covx_x = S*Covx_s;
end
if nargout >= 4
  % Compute entropy term:
  %
  % INT[ p(x_n, x_(n+1) | Y) log p(x_n | Y, x_(n+1)) dx_n dx_(n+1) ]
  Cov_joint = [Covx    Covx_x
               Covx_x' Covx_s];
  entropy = gaussian_entropy(logdet_cov(Cov_joint), size(Cov_joint,1)) ...
            - gaussian_entropy(logdet_cov(Covx_s), size(Covx_s,1));
                      
end
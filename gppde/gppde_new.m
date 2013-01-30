% [M_F, V_F] = GPPDE(X_Y, Y, S2_Y, THETA, X_G, G, ALPHA, D, X_F)

% Copyright (c) 2010 Jaakko Luttinen

% TODO: Each equation as a struct/cell: {X, y, K, alpha, D}
%
% X : input
% y : observations
% K : observation covariance matrix
% alpha : coefficients in (differential) equation
% D : differential operators in (differential) equation
%
% Then, noisy function observations are:
% { X_y, y, s2_y, 1, 0 }
% Noiseless derivative observations:
% { X_dy, dy, 1e-8, 1, 1 }
% E.g., heat equation:
% { X_g, zeros(..), 1e-8, [1; -a], [1 0; 0 2] }
%
% Actually, you should allow s2 to be a covariance matrix.  That way you
% could use other unknown functions in equations.. Or maybe even more
% generally, s2 could be a covariance function, thus its hyperparameters
% could be learnt simultaneously. :)
%
% And you can give several of these.
%
% And same for the predictions!

% TODO: Learn hyperparameters.

function [m_f, V_f] = gppde(varargin)

%
% Generate full joint covariance matrix for the equations
%


% Identity matrix for data noise
I_y = speye(size(X_y,1));
alpha0 = 1;
D0 = zeros(1, size(X_y,2));

% Variables for PDE numerical inaccuracy
I_g = speye(length(g));
s2_g = 1e-6;

% Joint data
yg = [y; g];

%
% DEBUG STUFF
%
% $$$ func = @(x) covfunc(x(:)', X_g(1,:), theta, alpha0,D0, alpha, D);
% $$$ mycheckgrad(func, 

% Joint prior covariance matrix for observations and PDE
K_g = covfunc(X_g, X_g, theta, alpha, D, alpha, D) + s2_g*I_g;
K_y = covfunc(X_y, X_y, theta, alpha0,D0, alpha0,D0) + s2_y*I_y;
K_y_g = covfunc(X_y, X_g, theta, alpha0,D0, alpha, D);
K_yg = [K_y, K_y_g; K_y_g', K_g];

% $$$ figure
% $$$ imagesc(K_yg)

L_yg = chol(K_yg, 'lower');

% Posterior covariance for predictive function values
K_y_f = covfunc(X_y, X_f, theta, alpha0,D0, alpha0,D0);
K_g_f = covfunc(X_g, X_f, theta, alpha,D, alpha0,D0);
K_yg_f = [K_y_f; K_g_f];
K_f = covfunc(X_f, X_f, theta, alpha0,D0, alpha0,D0);
Z = linsolve_tril(L_yg, K_yg_f);
V_f = K_f - Z' * Z;
m_f = Z' * linsolve_tril(L_yg, yg);


function K = covfunc(X1, X2, theta, alpha1, D1, alpha2, D2)

M = length(alpha1);
N = length(alpha2);
P = size(X1,2);

% For now, distances with unit length scales
D = bsxfun(@minus, ...
           reshape(X1, [size(X1,1),1,P]), ...
           reshape(X2, [1,size(X2,1),P]));

% Unit scale
K0 = theta(1)^2 * exp(-0.5/theta(2)^2*sum(D.^2,3));

K = 0;
for m=1:M
  for n=1:N
    % [c,E] = tmp(D1(m,:)+D2(n,:));
    
    % Evaluate exponents (E) and coefficients (c)
    b = D1(m,:)+D2(n,:) ;
    A = ntuples(floor(b/2));
    B = repmat(b, [size(A,1),1]);
    E = bsxfun(@minus, b, 2*A) ;
% $$$     c = (-1)^(sum(b)) * (-1).^(sum(A,2)) .* prod(npairsk(B, A), 2) % debug test
    c = (-1)^(sum(D1(m,:))) * (-1).^(sum(A,2)) .* prod(npairsk(B, A), 2) ;
    
    for l=1:length(c)
      Z = prod(bsxfun(@power, D, reshape(E(l,:), [1,1,P])), 3);
      Z = c(l) * Z;
      thetaexp = sum(0.5*(E(l,:)+b));
      K = K + alpha1(m)*alpha2(n) * theta(2)^(-2*thetaexp) * Z;
    end
  end
end
K = K.*K0;

% $$$ function [c, E] = term_coeff(d)
% $$$ 
% $$$ Z = ntuples(floor(d/2))
% $$$ D = repmat(d, [size(Z,1),1]);
% $$$ E = bsxfun(@minus, d, 2*Z);
% $$$ c = (-1).^(sum(Z,2)) .* prod(npairsk(D, Z), 2);


function test_gppde_1d

X_g = (-5:0.5:10)';
X_y = []; %(10*rand(1,2) - 5)';
X_f = (-5:0.1:10)';

% Data
theta = [10 1];
s2_y = theta(1)^2 * 1e-1;
y = []; %theta(1)*sin(2*pi*X_y(:)/theta(2)) + sqrt(s2_y)*randn(size(X_y,2),1);

% Define the partial differential equation
N_g = size(X_g, 1);
switch 2
 case 2   % f'' = -a*f
  D = [0; 2];
  alpha = [1; 0.5];
end
g = zeros(N_g,1);

% Inference: Get posterior predictive distribution
[m_f, V_f] = gppde(X_y, y, s2_y, theta, X_g, g, alpha, D, X_f);
figure
imagesc(V_f);

% Draw posterior predictive samples
I_f = speye(length(V_f));
s2_f = 1e-8; % numerical reasons
L_f = chol(V_f + s2_f*I_f, 'lower');
f = bsxfun(@plus, L_f * randn(length(L_f),10), m_f);
figure
plot(X_f, f);

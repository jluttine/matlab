
function test_gppde_3d

error('Not 3d yet')

% $$$ x_g = -5:0.5:5;
% $$$ [x1_g, x2_g] = meshgrid(
T = 10;
L = 10;
[x1_g, x2_g] = meshgrid(linspace(-0.5*T,1.5*T,40), linspace(-0.5*L,1.5*L,40));
X_g = [x1_g(:), x2_g(:)];
%X_g = [unifrnd(-0.5*T,1.5*T,[2000,1]), unifrnd(-0.5*L,1.5*L,[2000,1])];
X_y = [unifrnd(0,T, [0,1]), unifrnd(0,L, [0,1])];
[x1_f, x2_f] = meshgrid(linspace(0,T,20), linspace(0,L,20));
X_f = [x1_f(:), x2_f(:)];

% Gaussian process
theta = [10 1];

% Data
s2_y = theta(1)^2 * 1e-5;
y =  + sqrt(s2_y)*randn(size(X_y,1),1);

% Define the partial differential equation
N_g = size(X_g, 1);
D = [1 0
     0 2];
alpha = [1
         -1];
g = zeros(N_g,1);

% Inference: Get posterior predictive distribution
[m_f, V_f] = gppde(X_y, y, s2_y, theta, X_g, g, alpha, D, X_f);
% $$$ figure
% $$$ imagesc(V_f);

% Draw posterior predictive samples
samples = 9;
I_f = speye(length(V_f));
s2_f = 1e-8; % numerical reasons
L_f = chol(V_f + s2_f*I_f, 'lower');
f = bsxfun(@plus, L_f * randn(length(L_f),samples), m_f);
figure
for i = 1:samples
  subplot(ceil(sqrt(samples)), ceil(sqrt(samples)), i);
%  pcolor(reshape(f(:,i), size(x1_f)));
  contourf(x1_f, x2_f, reshape(f(:,i), size(x1_f)), 100);
  shading('flat')
end


function test_gppde_2d

% $$$ x_g = -5:0.5:5;
% $$$ [x1_g, x2_g] = meshgrid(
T = 1;
L = 1;
%[x1_g, x2_g] = meshgrid(linspace(-0*T,3*T,60), linspace(-1*L,2*L,40));
min_L = 0.0*L;
max_L = 1.0*L;
min_T = 0.0*T;
max_T = 1.0*T;
[x1_g, x2_g] = meshgrid(linspace(min_T,max_T,40), linspace(min_L,max_L,40));
X_g = [x1_g(:), x2_g(:)];
[x1_f, x2_f] = meshgrid(linspace(0,T,45), linspace(0,L,45));
X_f = [x1_f(:), x2_f(:)];



% Define the partial differential equation
N_g = size(X_g, 1);
switch 1
 case 1 
  %
  % HEAT EQUATION (d_t - d_x^2 = 0)
  %
  theta = [1 1/10]; % [scale, length scale]
  s2_y = theta(1)^2 * 1e-5;

  D = [1 0;     % differentiation
       0 2];
  alpha = [-1;   % coefficients
           0.02]; % (<- heat conduction)
  
  % Spatial boundary fixed to zero
  N_y = 30;
  bound1 = [linspace(min_T,max_T,N_y)', min_L*ones(N_y,1)];
  bound2 = [linspace(min_T,max_T,N_y)', max_L*ones(N_y,1)];
  X_y = [bound1; bound2];
  y = zeros(size(X_y,1),1);
  % Initial function at time T=min_T
  T_obs = min_T + 0.1*(max_T-min_T);
  bound_init = [T_obs*ones(N_y,1), linspace(min_L,max_L,N_y)'];
  X_y = [X_y; bound_init];
  y_init = 2*exp( -0.5*(bound_init(:,2)-0.5*(min_L+max_L)).^2 / 0.1^2);
  y = [y; y_init] + sqrt(s2_y)*randn(size(X_y,1),1);

  equation = 'Heat equation';
  
 case 2
  % 
  % WAVE EQUATION (d_t^2 - d_x^2 = 0)
  %
  theta = [1 1/30]; % [scale, length scale]
  s2_y = theta(1)^2 * 1e-4;

  D = [2 0;
       0 2];
  alpha = [-1;
           16]; % speed ^ 2
  % Spatial boundary fixed to zero
  N_y = 100;
  bound1 = [linspace(min_T,max_T,N_y)', min_L*ones(N_y,1)];
  bound2 = [linspace(min_T,max_T,N_y)', max_L*ones(N_y,1)];
  X_y = [bound1; bound2];
  y = zeros(size(X_y,1),1);
% $$$   % Initial function at time T=min_T
% $$$   bound_init = [min_T*ones(N_y,1), linspace(min_L,max_L,N_y)'];
% $$$   X_y = [X_y; bound_init];
% $$$   y_init = 2*exp( -0.5*(bound_init(:,2)-0.5*(min_L+max_L)).^2 / 0.03^2);
% $$$   y = [y; y_init];
  y = y + sqrt(s2_y)*randn(size(X_y,1),1);
  
  equation = 'Wave equation';

 case 3
  % 
  % LAPLACE EQUATION (d_x^2 + d_y^2 = 0)
  %
  theta = [1 1/20]; % [scale, length scale]
  s2_y = theta(1)^2 * 1e-5;

  D = [2 0;
       0 2];
  alpha = [1;
           1];
  X_y = zeros(0,2);
  y = zeros(0,1);
  
end
g = zeros(N_g,1);

% Inference: Get posterior predictive distribution
[m_f, V_f] = gppde(X_y, y, s2_y, theta, X_g, g, alpha, D, X_f);
% $$$ figure
% $$$ imagesc(V_f);

% Draw posterior predictive samples
samples = 1;
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
  map_colormap();
  cl = max( abs(f(:,i)) );
  set(gca, 'clim', [-cl, cl])
  title(equation)
  xlabel('time');
  xlabel('time');
  ylabel('location');
  set(gca, 'xtick', [], 'ytick', []);
end

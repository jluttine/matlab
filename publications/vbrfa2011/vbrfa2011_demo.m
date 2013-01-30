% VBRFA2011_DEMO  -  A simple comparison of four different noise models
%                    for VBPCA
%
% Noise models:
%  - Gaussian
%  - Multivariate Student-t
%  - Independent Student-t
%  - Laplace

% Last modified 2011-08-23
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function vbrfa2011_demo()

N = 50;
M = 2;
D = 1;

randn('state', 6);
rand('state', 6);

mu = [0;0];
W = [2; 1];
X = orth(randn(D,N)')' * sqrt(N);
s2 = 0.2;

% Generate normal data
Y = W*X + repmat(mu,1,N);

% Add noise
Yn = Y + sqrt(s2)*randn(M,N);

% Generate some outliers
Yno = Yn;
p = (rand(M,N) < 0.05);
Yno(p) = 5 + 8 * rand(sum(p(:)),1);

% Generate some missing values
Ynom = Yno;
pmv = (rand(M,N) < 0.0);
Ynom(pmv) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG: LAPLACE NOISE
% Ynom = Y + laplace_rand(0,1,M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Run different methods
%

Dh = D + 0; % number of components to estimate

for model=1:4

  % Isotropic noise
  noise_iso = noise_module_isotropic(M, 1, ...
                                     'init', struct('a_tau', 1e1, ...
                                                    'b_tau', 1e0), ...
                                     'update_tau', 1);
  
  W_module = factor_module_ard(Dh+1, M, ...
                               'prior', struct('a_alpha', 1e10, ...
                                               'b_alpha', 1e15), ...
                               'noise_module', noise_iso);
  
  X_module = factor_module_iid(Dh+1, N, ...
                               'prior', struct('mu', [1; zeros(Dh,1)], ...
                                               'CovX', diag([1e-6; ones(Dh,1)])));

  switch model
   case 1
    % Gaussian noise
    noise = noise_module_fixed(M, N, 1);
    method = 'gaussian';
    
   case 2

    % Multivariate t
    noise = noise_module_multivariate_t(M, N, ...
                                        'update_nu', 1, ...
                                        'init', struct('nu', 0.1));
    method = 'multi-t';

   case 3
    % Independent t
    noise = noise_module_independent_t(M, N, ...
                                       'update_nu', 1, ...
                                       'nu', 'pooled', ...
                                       'init', struct('nu', 0.1));
    method = 'ind-t';

   case 4
    % Laplace
    noise = noise_module_laplace(M, N);
    method = 'laplace';
    
  end

  Q = vbfa(D+1, Ynom, W_module, X_module, noise, ...
           'maxiter', 50, ...
           'rotate', true);
  Yh = Q.W'*Q.X;
  
  %
  % Show reconstructions
  %

  plot_results(Ynom, Yh)
  filename = ['/home/jluttine/papers/rpca_journal/figures/fig_demo_' method ...
              '.eps'];
  print('-depsc2', filename);

  rmse_recon(model) = rmse(Y-Yh);
  
end

rmse_recon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function plot_results(Y, Yh)
%subspace2d(Yh, [], Y)
figure
plot(Y(1,:), Y(2,:), 'x', 'markersize', 5, 'color', [.2 .2 .2]);
if ~isempty(Yh)
  hold on
  plot(Yh(1,:), Yh(2,:), 'ko', 'markersize', 3);
  nans = nan*ones(1,size(Y,2));
  plot([Y(1,:);Yh(1,:);nans], [Y(2,:);Yh(2,:);nans], ':', 'color', [.0 .0 .0]);
end
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'XTick', [], 'YTick', []);
set(gcf, 'units', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1), pos(2), 5 5])
set(gcf, 'PaperPositionMode', 'auto', 'paperunits', 'centimeters', 'PaperSize', [5 5]);

xmg = 0.05 * (max(Y(1,:)) - min(Y(1,:)));
ymg = 0.05 * (max(Y(2,:)) - min(Y(2,:)));
xl = [min(Y(1,:))-xmg, max(Y(1,:))+xmg];
yl = [min(Y(2,:))-ymg, max(Y(2,:))+ymg];
set(gca, 'xlim', xl, 'ylim', yl);

return


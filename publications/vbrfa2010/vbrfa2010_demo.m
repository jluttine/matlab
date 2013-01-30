% VBRFA2010_DEMO - A simple comparison of VB PCA and two robust extensions.
%
% The algorithms do not use ARD prior for the loadings W in order to keep
% the visualizations of the subspaces more clear. In addition, it helps
% avoiding pruning out relevant components too easily.
%
% NOTE: The robust algorithms might be less sensitive to the initialization
% if they estimate many components. With only one estimated component, they
% might kill all the components.
%
% This demo doesn't always show good results for robust algorithms. You may
% want to run a few times. In addition, sometimes there are several good
% interpretations of the data and the robust algorithms pick one.

% Last modified 2010-06-10
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function Q = vbrfa2010_demo

n = 50;
m = 2;
d = 1;
dh = d + 0; % number of components to estimate

randn('state', 6);
rand('state', 6);

mu = [0;0];
W = [2; 1];
X = orth(randn(d,n)')' * sqrt(n);
s2 = 0.2;

% Generate normal data
Y = W*X + repmat(mu,1,n);

% Add noise
Yn = Y + sqrt(s2)*randn(m,n);

% Generate some outliers
Yno = Yn;
p = (rand(m,n) < 0.05);
Yno(p) = 5 + 8 * rand(sum(p(:)),1);

% Generate some missing values
Ynom = Yno;
pmv = (rand(m,n) < 0.0);
Ynom(pmv) = NaN;

% Run different algorithms
options.init.tau = 1e2*ones(m,1); % the initialization of this can be crucial..
options.init.nu = 1;
options.prior.aw = 1e10*ones(d,1); % fix w by setting strong prior
options.prior.bw = 1e15*ones(d,1); % fix w by setting strong prior
options.update_nu = 1;
options.update_w = 1;
options.rotate = 0;
options.common_nu = false;
options.common_tau = true
options.maxiter = 50;

% PCA
disp('Run PCA')
options.robustness = 'none';
results_p = vbrfa(Ynom, dh, options);
W_p = results_p.W;
X_p = results_p.X;
S_p = results_p.Sv;
mu_p = results_p.mu;
nu_p = results_p.nu;

% Robust PCA with multivariate Student-t
disp('Run robust PCA with multivariate Student-t')
options.robustness = 'multivariate-t';
results_rp = vbrfa(Ynom, dh, options);
W_rp = results_rp.W;
X_rp = results_rp.X;
S_rp = results_rp.Sv;
mu_rp = results_rp.mu;
nu_rp = results_rp.nu;

% Robust PCA with independent Student-t
disp('Run robust PCA with independent Student-t')
options.robustness = 'independent-t';
results_rvb = vbrfa(Ynom, dh, options);
W_rvb = results_rvb.W;
X_rvb = results_rvb.X;
S_rvb = results_rvb.Sv;
mu_rvb = results_rvb.mu;
nu_rvb = results_rvb.nu;

% Reconstruct
Y_p = W_p*X_p + repmat(mu_p,1,n);
Y_rp = W_rp*X_rp + repmat(mu_rp,1,n);
Y_rvb = W_rvb*X_rvb + repmat(mu_rvb,1,n);

plot_results(Ynom, []);
%print('-depsc2', '/home/jluttine/papers/vbrfa/figures_journal/fig_demo_data.eps');
plot_results(Ynom, Y_p)
%print('-depsc2', '/home/jluttine/papers/vbrfa/figures_journal/fig_demo_pca.eps');
plot_results(Ynom, Y_rp)
%print('-depsc2', '/home/jluttine/papers/vbrfa/figures_journal/fig_demo_multi-t.eps');
plot_results(Ynom, Y_rvb)
%print('-depsc2', '/home/jluttine/papers/vbrfa/figures_journal/fig_demo_indep-t.eps');

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

% $$$ % Mark outliers
% $$$ hold on
% $$$ indeces = sum(p,1)>0;
% $$$ plot(VYno(1,indeces), VYno(2, indeces), 'go', 'MarkerSize', 8);

% $$$ % Mark observations with missing values
% $$$ indeces = sum(pmv,1)>0;
% $$$ scatter(VYno(1,indeces), VYno(2, indeces), 'yo');

return


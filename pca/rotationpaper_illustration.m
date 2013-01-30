function rotationpaper_illustration

% (prior) variances of y,k,x
vy = 0.2^2
vk = 1^2
vx = 1^2

% observation
y = 0.8;

% posterior parameters
v_k = 0.4;
m_k = 0.2;
v_x = 0.4;
m_x = 0.2;

figsize = 6; % centimeters

% grid approximate the true posterior
[K,X] = meshgrid(-2.5:0.01:2.5, -2.5:0.01:2.5);
logPykx = norm_lpdf(y,K.*X,sqrt(vy)) + norm_lpdf(K,0,sqrt(vk)) + ...
       norm_lpdf(X,0,sqrt(vx));
%Py = sum(exp(logPykx(:)));

% $$$ %levels = [0.3;0.2;0.1];
% $$$ figure
% $$$ for ind=1:levels
% $$$   contour(K,X,exp(logPykx));
% $$$   hold on
% $$$ end

levels = [0.15 0.45 0.75];
figure
Pykx = exp(logPykx);
Pykx = Pykx / max(Pykx(:));
contour(K,X,Pykx, levels);
xl = xlim;
yl = ylim;
hold on

muQ = [m_x;m_k];
CovQ = diag([v_x;v_k]);
plot_product_frame(K,X,Pykx);
set(gca, 'xlim', xl, 'ylim', yl);
% $$$ set(gca, 'xtick', [], 'ytick', []);
xlabel('$x$')
ylabel('$w$')
adjust_figure
print2latex('/home/jluttine/clim/rotateVBPCA/fig_banana_posterior', figsize);

muQ = [m_x;m_k];
CovQ = diag([v_x;v_k]);
plot_product_frame(K,X,Pykx,muQ,CovQ);
set(gca, 'xlim', xl, 'ylim', yl);
% $$$ set(gca, 'xtick', [], 'ytick', []);
xlabel('$x$')
ylabel('$w$')
adjust_figure
print2latex('/home/jluttine/clim/rotateVBPCA/fig_banana_init', figsize);

%return


% $$$ % Analytic posterior (for sum model)
% $$$ w = ones(2,1); % corresponds to summing x and k
% $$$ Cov = inv(1/vy*w*w' + diag([vx,vk].^(-1)));
% $$$ mu = Cov * 1/vy * w * y;
% $$$ 
% $$$ muQ = [m_x;m_k];
% $$$ CovQ = diag([v_x;v_k]);
% $$$ plot_frame(mu,Cov,muQ,CovQ);
% $$$ set(gca, 'xtick', [], 'ytick', []);
% $$$ xl = xlim;
% $$$ yl = ylim;
% $$$ xlabel('$x_1$')
% $$$ ylabel('$x_2$')

for ind=1:100
  % Update k
  [m_k, v_k] = update(y,vy, vk, m_x,v_x);
  if ind <=16
    muQ_old = muQ;
    CovQ_old = CovQ;
    muQ = [m_x;m_k];
    CovQ = diag([v_x;v_k]);
    %plot_frame(mu,Cov,muQ,CovQ,muQ_old,CovQ_old);
    plot_product_frame(K,X,Pykx,muQ,CovQ,muQ_old,CovQ_old);
    set(gca, 'xlim', xl, 'ylim', yl);
% $$$     set(gca, 'xtick', [], 'ytick', []);
    xlabel('$x$')
    ylabel('$w$')
    adjust_figure
    print2latex(sprintf('/home/jluttine/clim/rotateVBPCA/fig_banana_iter%da', ...
                        ind), figsize);
  end
  % Update x
  [m_x, v_x] = update(y,vy, vx, m_k,v_k);
  if ind <=16
    muQ_old = muQ;
    CovQ_old = CovQ;
    muQ = [m_x;m_k];
    CovQ = diag([v_x;v_k]);
%    plot_frame(mu,Cov,muQ,CovQ,muQ_old,CovQ_old);
    plot_product_frame(K,X,Pykx,muQ,CovQ,muQ_old,CovQ_old);
    set(gca, 'xlim', xl, 'ylim', yl);
% $$$     set(gca, 'xtick', [], 'ytick', []);
    xlabel('$x$')
    ylabel('$w$')
    adjust_figure
    print2latex(sprintf('/home/jluttine/clim/rotateVBPCA/fig_banana_iter%db', ...
                        ind), figsize);
  end
  fprintf('yh = %.3f+/-%.3f\n', m_x*m_k, nan)
% $$$   fprintf('yh = %.3f+/-%.3f\n', m_x+m_k, sqrt(v_x+v_k))
end

muQ = [m_x;m_k];
CovQ = diag([v_x;v_k]);
%plot_frame(mu,Cov,muQ,CovQ);
plot_product_frame(K,X,Pykx,muQ,CovQ);
set(gca, 'xlim', xl, 'ylim', yl);
% $$$ set(gca, 'xtick', [], 'ytick', []);
xlabel('$x$')
ylabel('$w$')
adjust_figure
print2latex('/home/jluttine/clim/rotateVBPCA/fig_banana_converge', figsize);

m_x
m_k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_frame(mu, Cov, muQ, CovQ,muQ_old,CovQ_old)
figure
draw_gaussian(mu, Cov, 1:3, 'k', 'linewidth', 3);
hold on
draw_gaussian(muQ, CovQ, 1:3, 'r', 'linewidth', 2);
if nargin == 6
  draw_gaussian(muQ_old, CovQ_old, 1:3, '--', 'color', [1 0.4 0.4], ...
                'linewidth', 1);
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_product_frame(K,X,Pykx,muQ,CovQ,muQ_old,CovQ_old)
levels = 3 %[0.15 0.45 0.75];
%levels = 4;

figure
contour(K,X,Pykx, levels, 'k', 'linewidth', 1);

if nargin >= 4
  Q = mvnpdf([K(:),X(:)], muQ', CovQ);
  Q = reshape(Q / max(Q(:)), size(K));
  hold on
  if nargin == 7
    oldQ = mvnpdf([K(:),X(:)], muQ_old', CovQ_old);
    oldQ = reshape(oldQ / max(oldQ(:)), size(K));
    hold on
    contour(K,X,oldQ, levels, '--', 'linecolor', [0.9 0.65 0.65], 'linewidth', 1);
  end
  contour(K,X,Q, levels, 'r-', 'linewidth', 1);
% $$$ draw_gaussian(muQ, CovQ, 1:3, 'r', 'linewidth', 2);
% $$$ if nargin == 7
% $$$   draw_gaussian(muQ_old, CovQ_old, 1:3, '--', 'color', [1 0.4 0.4], ...
% $$$                 'linewidth', 1);
% $$$ end  
end

function [m_k, v_k] = update(y, vy, vk, m_x, v_x)

% Sum model
% $$$ v_k = 1 / ( 1/vy + 1/vk );
% $$$ m_k = v_k * 1/vy * (y-m_x);

% Product model
v_k = 1 / ( (m_x^2+v_x)/vy + 1/vk );
m_k = v_k * m_x/vy * y;

function adjust_figure
len = 2;
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'TickDir', 'out');
set(gca, 'Box', 'off' );
% $$$ set(gca, 'Units', 'centimeters');
% $$$ pos = get(gca, 'Position');
% $$$ set(gca, 'Position', [pos(1:2), 8 8]);
set(gca, 'xtick', -2:2);
set(gca, 'ytick', -2:2);
% $$$ set(gcf, 'PaperUnits', 'centimeters');
% $$$ set(gcf, 'Units', 'centimeters');
% $$$ set(gcf, 'PaperSize', [len len]);
% $$$ pos = get(gcf, 'PaperPosition');
% $$$ set(gcf, 'PaperPosition', [pos(1:2) len len]);
% $$$ pos = get(gcf, 'Position');
% $$$ set(gcf, 'Position', [pos(1:2) len len]);

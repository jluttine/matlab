function aistats2012_slides_gp

randn('state', 1);

plot_gp(true, 0.2, 0.4, 100, false, ...
        '/home/jluttine/papers/aistats2012/slides/fig_short', ...
        8, 8);
plot_gp(true, 1.5, 0.4, 100, false, ...
        '/home/jluttine/papers/aistats2012/slides/fig_long', ...
        8, 8);
plot_gp(false, 1, 0.2, 6, 5, ...
        '/home/jluttine/papers/aistats2012/slides/fig_gp', ...
        6, 4);

function plot_gp(noisy, lengthscale, noise, N_all, samples, filename, width, height);
xmin = 0;
xmax = 10;

Nh = 1000;
xh = linspace(xmin, xmax, Nh);

x_all = linspace(xmin+0.5, xmax-0.5, N_all);
%x_all = unifrnd(xmin, xmax, 1, N_all);

% Generate data
covfunc_all = gp_cov_sum(gp_cov_se(sq_dist(x_all), ...
                                 'lengthscale', lengthscale), ...
                       gp_cov_scale(gp_cov_delta(N_all), ...
                                    'scale', noise));

K_y = covfunc_all([]);
L_y = chol(K_y + 1e-6*eye(N_all), 'lower');
y_all = L_y * randn(N_all,1);

for N = [N_all]
  x = x_all(:,1:N);
  y = y_all(1:N);

  covfunc_y = gp_cov_sum(gp_cov_se(sq_dist(x), ...
                                   'lengthscale', lengthscale), ...
                         gp_cov_scale(gp_cov_delta(N), ...
                                      'scale', noise));
  covfunc_yf = gp_cov_se(sq_dist(x,xh), ...
                         'lengthscale', lengthscale);
  if noisy
    covfunc_f = gp_cov_sum(gp_cov_se(zeros(Nh,1), ...
                                     'lengthscale', lengthscale), ...
                           gp_cov_scale(gp_cov_wrap(ones(Nh,1)), ...
                                        'scale', 0.3));
  else
    covfunc_f = gp_cov_se(zeros(Nh,1), ...
                          'lengthscale', lengthscale);
  end

  K_y = covfunc_y([]);
  K_yf = covfunc_yf([]);
  k_f = covfunc_f([]);

  [fh_mean, fh_var] = gp_predict(y, K_y, K_yf, k_f);
  figure
  errorplot(xh, fh_mean, sqrt(fh_var))
  hold on
  plot(x, y, 'r+')
  set(gca, 'xtick', [], 'ytick', []);
  set(gca, 'box', 'off');
  set_figure_size(width,height);
  print('-dpdf', [filename, '.pdf'])
  yl = get(gca, 'ylim');

  if samples
    covfunc_f = gp_cov_se(sq_dist(xh), ...
                          'lengthscale', lengthscale);
    K_f = covfunc_f([]);
    K_post = K_f - K_yf'*(K_y\K_yf) + 1e-6*eye(Nh);
    L_post = chol(K_post, 'lower');
    f = bsxfun(@plus, ...
               K_yf'*(K_y\y), ...
               L_post*randn(Nh,samples));
    figure
    plot(xh, f);
    hold on
    plot(x, y, 'k+');
    set(gca, 'xtick', [], 'ytick', []);
    set(gca, 'box', 'off');
    set(gca, 'ylim', yl);
    set_figure_size(width,height);  
    print('-dpdf', [filename, '_samples.pdf'])
  end
   
end


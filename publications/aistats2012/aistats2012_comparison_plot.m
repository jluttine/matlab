function aistats2012_plot_comparison()

% Results for F
N_chol = [10^3.0 10^3.5 10^4.0 10^4.5 10^5.0  10^5.5];
t_chol = [  2.53  30.05 362.22 2859.2 27451.8  NaN  ];
% $$$ t_chol = [  0.96   10.4    151   1078  11036];
N_kron = [10^3.0 10^3.5 10^4.0 10^4.5 10^5.0 10^5.5 10^6.0 10^7.0 10^8.0];
t_kron = [0.026   0.046  0.137  0.565   2.01   7.13  20.15  358.8   4017];
% $$$ t_kron = [0.012    0.06   0.56   1.09   1.93   15.2    108   1200];

figure
loglog(N_chol, t_chol, 'r+-');
hold on
loglog(N_kron, t_kron*1e2, 'bo-');
loglog(N_kron, t_kron, 'bo--');
hl = legend({'Cholesky', 'Kronecker (100)', 'Kronecker (1)'}, ...
       'Location', 'SouthEast');
% $$$ hl = legend({'Cholesky (exact)', 'Kronecker (100 samples)', 'Kronecker (1 sample)'}, ...
% $$$        'Location', 'SouthEast');
hx = xlabel('Number of data points');
hy = ylabel('CPU time (seconds)');
set(gca, 'FontSize', 10);
set(hl, 'FontSize', 8);
set(hx, 'FontSize', 10);
set(hy, 'FontSize', 10);
set(gca, 'xlim', [1e3 1e8])
set(gca, 'xtick', [1e3 1e4 1e5 1e6 1e7 1e8])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'off')
set(gca, 'ylim', [1e-2, 1e6])
decades_equal(gca);
%set(gca, 'Box', 'off' );
set_figure_size(10,11)
print('-depsc2', '/home/jluttine/papers/gpkron/figures/fig_artificial_f');

% Results for theta
N_chol = [10^(3.0) 10^(3.5) 10^(4.0) 10^(4.5) 10^(5.0)];
t_chol = [  25.221  100.168  572.920 2490.262     NaN];
t_chol2 = [  8.247   40.374  170.632  873.243     NaN];
% $$$ t_chol = [   7.36     36.5     96.9      956     4092];
% $$$ t_chol2 = [  4.90     19.0     69.4      333     1605];
N_kron = [ 10^3.0 10^3.5 10^4.0 10^4.5 10^5.0  10^6.0   10^7.0 10^8.0];
t_kron = [  6.274 11.739 64.465 60.628 241.449 3365.156  63264    NaN];
t_kron2 =[  0.047  0.131  0.238  0.785   2.443   27.045    450    NaN];
% $$$ t_kron = [    0.9    5.6   37.4   50.7    125   1017   5801    NaN];
% $$$ t_kron2 =[   0.03   0.13    0.8    1.5    2.5   18.5    133   1500];

figure
loglog(N_chol, t_chol, 'r+-');
hold on
loglog(N_chol, t_chol2, 'r+--');
loglog(N_kron, t_kron, 'bo-');
loglog(N_kron, t_kron2, 'bo--');
hl = legend({'Cholesky (eff.)', 'Cholesky (single)', ...
             'Kronecker (eff.)',  'Kronecker (single)'}, ...
       'Location', 'SouthEast');
hx = xlabel('Number of data points');
hy = ylabel('CPU time (seconds)');
set(gca, 'FontSize', 10);
set(hl, 'FontSize', 8);
set(hx, 'FontSize', 10);
set(hy, 'FontSize', 10);
set(gca, 'xlim', [1e3 1e7])
set(gca, 'ylim', [1e-2, 1e5])
set(gca, 'xtick', [1e3 1e4 1e5 1e6 1e7])
set(gca, 'xminortick', 'off')
set(gca, 'yminortick', 'off')
decades_equal(gca);
%set(gca, 'Box', 'off' );
set_figure_size(10,11)
print('-depsc2', '/home/jluttine/papers/gpkron/figures/fig_artificial_theta');

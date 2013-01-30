function nips2011_plot_comparison()

% Results for F
N_chol = [10^3.0 10^3.5 10^4.0 10^4.5 10^5.0];
t_chol = [  0.96   10.4    151   1078  11036];
%t_chol = [   1.8   14.4    166   1187   5506];
% $$$ N_chol = [1e3 1e4 10^(4.3) 10^(4.5) 10^(4.8)];
% $$$ t_chol = [1.8 173     554      1238     3290];
N_kron = [10^3.0 10^3.5 10^4.0 10^4.5 10^5.0 10^6.0 10^7.0 10^8.0];
t_kron = [0.012    0.06   0.56   1.09   1.93   15.2    108   1200];

figure
loglog(N_chol, t_chol, 'r+-');
hold on
loglog(N_kron, t_kron*1e2, 'bo-');
loglog(N_kron, t_kron, 'bo--');
hl = legend({'Cholesky (exact)', 'Kronecker (100 samples)', 'Kronecker (1 sample)'}, ...
       'Location', 'SouthEast');
hx = xlabel('Number of data');
hy = ylabel('CPU time (seconds)');
%title('Time for inference of F');
%set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'FontSize', 12);
set(hl, 'FontSize', 11);
set(hx, 'FontSize', 12);
set(hy, 'FontSize', 12);
set(gca, 'xlim', [1e3 1e8])
%set(gca, 'ylim', [10^(-2) 1e5])
set(gca, 'Box', 'off' );
set_figure_size(10,10)
print('-depsc2', '/home/jluttine/papers/gpkron/fig_artificial_f');

% Results for theta
N_chol = [10^(3.0) 10^(3.5) 10^(4.0) 10^(4.5) 10^(5.0)];
t_chol = [   7.36     36.5     96.9      956     4092];
t_chol2 = [  4.90     19.0     69.4      333     1605];
N_kron = [ 10^3.0 10^3.5 10^4.0 10^4.5 10^5.0 10^6.0 10^7.0 10^8.0];
t_kron = [    0.9    5.6   37.4   50.7    125   1017   5801    NaN];
t_kron2 =[   0.03   0.13    0.8    1.5    2.5   18.5    133   1500];

figure
loglog(N_chol, t_chol, 'r+-');
hold on
loglog(N_chol, t_chol2, 'r+--');
loglog(N_kron, t_kron, 'bo-');
loglog(N_kron, t_kron2, 'bo--');
hl = legend({'Cholesky (effective)', 'Cholesky (single)', ...
             'Kronecker (effective)',  'Kronecker (single)'}, ...
       'Location', 'SouthEast');
hx = xlabel('Number of data');
hy = ylabel('CPU time (seconds)');
%title('Sampling time of \theta')
%set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'FontSize', 12);
set(hl, 'FontSize', 11);
set(hx, 'FontSize', 12);
set(hy, 'FontSize', 12);
set(gca, 'xlim', [1e3 1e7])
set(gca, 'Box', 'off' );
set_figure_size(10,10)
print('-depsc2', '/home/jluttine/papers/gpkron/fig_artificial_theta');

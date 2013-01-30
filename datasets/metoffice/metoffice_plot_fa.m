function metoffice_plot_fa(Q)

[D,N] = size(Q.X);
M = size(Q.W,2);

E2 = reshape(Q.CovX, [D^2, N]);
E = sqrt(E2(1:(D+1):end,:));

hax = tserrorplot(Q.X',2*E');
ylim_centered(hax);
ylim_scale(hax, 1.2);

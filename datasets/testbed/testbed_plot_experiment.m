function testbed_plot_experiment(Q)
% testbed_plot_experiment(Q)
%
% Plots latent signals, spatial components and reconstructions.

% Plot latent signals
testbed_states(Q);

% Plot spatial components
testbed_loadings(Q, 'resolution', [50 30], 'method', 'rbf');

% Plot reconstructions of problematic stations
testbed_reconstructions(Q, 'stations', [12 13 20 21 28 48 57 59]);

% $$$ M = size(Q.W,1);
% $$$ group = 10;
% $$$ for m=1:group:M
% $$$   stations = m:min(m+group-1, M);
% $$$   testbed_rcplot(Q, 'stations', stations);
% $$$ end

% VBRFA2011_TESTBED_PLOT - Make the Testbed plots for the VBRFA journal
% paper.

% Last modified 2011-06-15
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function data = vbrfa2011_testbed_plot(model)

% Result files
folder = '/share/climate/jluttine/testbed';
%files = {'testbed_results_ind-t_D=30_20110616'};
files = {'testbed_results_gaussian_D=30_20110617', ...
         'testbed_results_multi-t_D=30_20110617', ...
         'testbed_results_ind-t_D=30_20110617', ...
         'testbed_results_laplace_D=30_20110617'};

% Pick stations from: 12 13 20 21 28 48 57 59
corrupted = [12 13 20 21 28 48 57 59];
stations = [5 corrupted([8 7 6 1 3])];

%
% PLOT DATA OBSERVATIONS
%

data = testbed_loaddata();
data = testbed_preprocess(data);
return
plot_timeseries(data.time, ...
                data.observations(stations,:));
%return

print('-depsc2', '-loose', ...
      '/home/jluttine/papers/rpca_journal/fig_testbed_timeseries_data');
%
% PLOT RESULTS FOR EACH METHOD
%

suffix = {'gaussian', 'multi-t', 'ind', 'laplace'};

% Plot
if nargin < 1
  model = 1:length(files);
end
for n=model
  fprintf('Loading results from file:\n%s\n', files{n});
  Q = load(files{n}, 'W', 'X', 'W_struct');
% $$$   plot_loadings(Q,suffix{n});
  plot_reconstructions(Q,stations,suffix{n});
  %figure, semilogy(Q.W_struct.alpha)
  %fprintf('alpha(%d): %.4e\n', [1:length(Q.w); Q.w(:)']);
end

% $$$ %
% $$$ % PLOT THE STATIONS
% $$$ %
% $$$ 
% $$$ fig = figure();
% $$$ 
% $$$ % Plot the research area on the map of Finland
% $$$ hax(1) = subplot(1,2,1);
% $$$ testbed_plot_finland();
% $$$ 
% $$$ % Plot the stations on a topographic map
% $$$ hax(2) = subplot(1,2,2);
% $$$ set(hax(2), 'FontSize', 9);
% $$$ testbed_plot_topography();
% $$$ hold on;
% $$$ testbed_coast();
% $$$ map_grid();
% $$$ testbed_plot_stations(data.coordinates);
% $$$ hcb = colorbar('peer',hax(2),'FontSize',9);
% $$$ 
% $$$ % Set layout
% $$$ bottom = 0.07;
% $$$ top = 0.92;
% $$$ set(hax(1), 'Units','Normalized', 'Position',[0.0 bottom 0.22 (top-bottom)]);
% $$$ set(hax(2), 'Units','Normalized', 'Position',[0.27 bottom 0.63 (top-bottom)]);
% $$$ set(hcb, 'Units','Normalized', 'Position',[0.91 bottom 0.02 (top-bottom)]);
% $$$ set_figure_size(15,6.3,fig);
% $$$ 
% $$$ % $$$ % Print
% $$$ % $$$ print('-depsc2', '-loose', ...
% $$$ % $$$       '/home/jluttine/papers/vbrfa/figures_journal/fig_testbed_stations');
% $$$     
% $$$ 

function plot_timeseries(time,Y)

% Plot timeseries
hax = tsplot(time, Y, 'k');

% Set the range of the x-axis
N = length(hax);
set(hax, 'xlim', [min(time), max(time)]);

% Put date ticks on x-axis
dateticks(hax, 'x', 'mmmyy', 'keeplimits');
set(hax(1:(N-1)), 'xticklabel', []);
set_ticks_fontsize(hax,9);
xtick = get(hax(end), 'xtick');
set(hax, 'xtick', xtick(2:end));

% Remove ticks from y-axis
set(hax, 'ylim', [-50 50], 'yticklabel', [], 'ytick', []);

% Set the layout of the subplots
set_subplot_positions(hax, N, 1, [0.01 0.01 0.01 0.07], [0.02 0.02]);
set_figure_size(15,9)

function plot_reconstructions(Q,stations,suffix)

% Reconstructions
Yh = Q.W(:,stations)' * Q.X;

% Remove time instances with no observations to make the plot more clear
data = load(['/share/climate/jluttine/testbed/' ...
             'testbed_vbrfa2011_traindata']);
% $$$ data = testbed_loaddata();
% $$$ data = testbed_preprocess(data);
ind = sum(~isnan(data.observations),1) == 0;
Yh(:,ind) = nan;

% Plotplotplot
plot_timeseries(data.time, Yh);

% Print figure
print('-depsc2', '-loose', ...
      ['/home/jluttine/papers/rpca_journal/fig_testbed_timeseries_', ...
       suffix]);

function plot_loadings(Q,suffix)

%
% Plot spatial components
%

data = testbed_loaddata();
data = testbed_preprocess(data);

components = 1:4;
hax = testbed_loadings(Q.W(components,:)', ...
                       data.coordinates, ...
                       'resolution', [200 150], ...
                       'method', 'rbf');
% Subplot layout
set_subplot_positions(hax, 2, 2, [0.01 0.01 0.08 0.06], [0.08 0.1]);
% Colorbar layout
colorbars(hax, ...
          'Size', 0.02, ...
          'Separation', 0.0, ...
          'Location', 'EastOutside', ...
          'FontSize', 9);
% X-labels (and layout)
for n=1:numel(components)
  xlabels(hax(n),sprintf('(%d)', components(n)), 'Separation',0.03, 'FontSize', 9);
end
% Print figure
set_figure_size(15,8.5)
% $$$ print('-depsc2', '-loose', ...
% $$$       ['/home/jluttine/papers/vbrfa/figures_journal/fig_testbed_loadings_', ...
% $$$        suffix]);

% Plot reconstructions of problematic stations
%testbed_rcplot(Q, 'stations', [12 13 20 21 28 48 57 59]);


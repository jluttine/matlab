function nips2011_plot_mohsst5()

dir = '/home/jluttine/matlab/publications/nips2011/';

GP1 = load([dir 'results_nips2011_mohsst5_gpkron_uniform_20110528']);
GP2 = load([dir 'results_nips2011_mohsst5_gpkron_pattern_20110528']);
PCA1 = load([dir 'results_nips2011_mohsst5_vbpca_uniform_20110601']);
PCA2 = load([dir 'results_nips2011_mohsst5_vbpca_pattern_20110601']);

Y = PCA2.Ytest;
I = ~isnan(PCA2.Ytrain);
Y(I) = PCA2.Ytrain(I);

data = mohsst5_loaddata();
[LON,LAT] = meshgrid(data.longitude, data.latitude);
weights = 1./cosd(LAT(:));
weights = mohsst5_remove_land(weights);

%
% Means (uniform)
%

t = 800;
% t = 1200;

cl = max(max(abs(Y(:,t))), ...
         max(max(abs(GP1.F(:,t))), ...
             max(abs(PCA1.F(:,t)))));

dir = '/home/jluttine/papers/gpkron/';

figure
mohsst5_mapplot(GP1.F(:,t));
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_gp_mean']);


figure
mohsst5_mapplot(PCA1.F(:,t))
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_pca_mean']);

%
% DATA
%

figure
mohsst5_mapplot(Y(:,t))
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_data']);

%
% STD (uniform)
%

GPstd = sqrt(GP1.stdF(:,t).^2 + weights*mean(GP1.results.theta(4,1000:end).^2));
PCAstd = sqrt(PCA1.stdF(:,t).^2 + 1./PCA1.results.Tau(:,t));

cl = max( max(GPstd), max(PCAstd) );

figure
mohsst5_mapplot(GPstd);
colormap_scale();
polish(gca);
set(gca, 'clim', [0 cl]);
print('-depsc2', [dir 'fig_sst_gp_std']);

figure
mohsst5_mapplot(PCAstd);
colormap_scale();
polish(gca);
set(gca, 'clim', [0 cl]);
print('-depsc2', [dir 'fig_sst_pca_std']);

disp(datestr(data.time(t)));

%return

%return

%
% Means (pattern)
%

t = 1429;

cl = max( max(abs(GP2.F(:,t))), max(abs(PCA2.F(:,t))) );

% GP

figure
mohsst5_mapplot(GP2.F(:,t));
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_pattern_gp_mean']);

% PCA

figure
mohsst5_mapplot(PCA2.F(:,t))
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_pattern_pca_mean']);

% Data

figure
mohsst5_mapplot(Y(:,t))
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_pattern_data']);

figure
mohsst5_mapplot(PCA2.Ytrain(:,t))
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
print('-depsc2', [dir 'fig_sst_pattern_train']);

disp(datestr(data.time(t)));

% $$$ figure
% $$$ EGP = abs(Y-GP2.F);
% $$$ EPCA = abs(Y-PCA2.F);
% $$$ mohsst5_mapplot(EGP(:,t)-EPCA(:,t))
% $$$ colormap_redblue();
% $$$ polish(gca);
% $$$ %set(gca, 'clim', [-cl cl]);

% $$$ figure
% $$$ mohsst5_mapplot(E2(:,t))
% $$$ colormap_scale();
% $$$ polish(gca);
%set(gca, 'clim', [-cl cl]);

%
% STD (pattern)
%

GPstd = sqrt(GP2.stdF(:,t).^2 + weights*mean(GP2.results.theta(4,1000:end).^2));
PCAstd = sqrt(PCA2.stdF(:,t).^2 + 1./PCA2.results.Tau(:,t));

cl = max( max(GPstd), max(PCAstd) );
%cl = max( max(GP2.stdF(:,t)), max(PCA2.stdF(:,t)) );

figure
mohsst5_mapplot(GPstd);
colormap_scale();
polish(gca);
set(gca, 'clim', [0 cl]);
print('-depsc2', [dir 'fig_sst_pattern_gp_std']);

figure
mohsst5_mapplot(PCAstd);
colormap_scale();
polish(gca);
set(gca, 'clim', [0 cl]);
print('-depsc2', [dir 'fig_sst_pattern_pca_std']);

return

function polish(hax)
hcb = colorbar('SouthOutside');
set(hcb, 'Units', 'normalized');
set_colorbar_position(hcb, hax, 0.03, 0.02);
set_figure_size(10,8);

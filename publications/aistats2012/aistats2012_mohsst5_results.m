function aistats2012_mohsst5_results()


%
% DATA
%

disp('--- Results for UNIFORM case ---')
show_results('data', 1);
disp('Kronecker GP')
show_results('gpkron', 1, 'results_aistats2012_mohsst5_gpkron_uniform_20110528');
disp('Kronecker GP (2 covfuncs)')
show_results('gpkron2', 1, 'results_aistats2012_mohsst5_gpkron-cov2_uniform_20111010');
disp('VB PCA')
show_results('pca', 1, 'results_aistats2012_mohsst5_vbpca_uniform_20111007');
disp('GPFA')
show_results('gpfa', 1, 'results_aistats2012_mohsst5_gpfa_uniform_20111006');

disp('--- Results for PATTERN case ---')
show_results('data', 2);
disp('Kronecker GP')
show_results('gpkron', 2, 'results_aistats2012_mohsst5_gpkron_pattern_20110528');
disp('Kronecker GP (2 covfuncs)')
show_results('gpkron2', 2, 'results_aistats2012_mohsst5_gpkron-cov2_pattern_20111010');
disp('VB PCA')
show_results('pca', 2, 'results_aistats2012_mohsst5_vbpca_pattern_20111007');
disp('GPFA')
show_results('gpfa', 2, 'results_aistats2012_mohsst5_gpfa_pattern_20111006');

return

% $$$ GP1 = load([dir 'results_aistats2012_mohsst5_gpkron-cov2_uniform_20110930']);
% $$$ % $$$ GP1 = load([dir 'results_aistats2012_mohsst5_gpkron-cov2_uniform_20110907']);
% $$$ % $$$ GP2 = load([dir 'results_aistats2012_mohsst5_gpkron-cov2_pattern_20110908']);
% $$$ PCA1 = load([dir 'results_aistats2012_mohsst5_vbpca_uniform_20110601']);
% $$$ PCA2 = load([dir 'results_aistats2012_mohsst5_vbpca_pattern_20110601']);
% $$$ GPFA1 = load([dir 'results_aistats2012_mohsst5_gpfa_uniform_20110921']);
% $$$ GPFA2 = load([dir 'results_aistats2012_mohsst5_gpfa_pattern_20111001']);
% $$$ 
% $$$ [gp1_rmse_train, gp1_rmse_test] = aistats2012_mohsst5_rmse(GP1.F, 1);
% $$$ fprintf('GP, uniform: train=%.4f and test=%.4f\n', gp1_rmse_train, gp1_rmse_test);
% $$$ [pca1_rmse_train, pca1_rmse_test] = aistats2012_mohsst5_rmse(PCA1.F, 1);
% $$$ fprintf('PCA, uniform: train=%.4f and test=%.4f\n', pca1_rmse_train, pca1_rmse_test);
% $$$ [pca2_rmse_train, pca2_rmse_test] = aistats2012_mohsst5_rmse(PCA2.F, 2);
% $$$ fprintf('PCA, pattern: train=%.4f and test=%.4f\n', pca2_rmse_train, pca2_rmse_test);
% $$$ [gpfa1_rmse_train, gpfa1_rmse_test] = aistats2012_mohsst5_rmse(GPFA1.F, 1);
% $$$ fprintf('GPFA, pattern: train=%.4f and test=%.4f\n', gpfa2_rmse_train, gpfa2_rmse_test);
% $$$ [gpfa2_rmse_train, gpfa2_rmse_test] = aistats2012_mohsst5_rmse(GPFA2.F, 2);
% $$$ fprintf('GPFA, pattern: train=%.4f and test=%.4f\n', gpfa2_rmse_train, gpfa2_rmse_test);
return

% $$$ Y = GP1.Ytest;
% $$$ I = ~isnan(GP1.Ytrain);
% $$$ Y(I) = GP1.Ytrain(I);
% $$$ 
% $$$ data = mohsst5_loaddata();
% $$$ [LON,LAT] = meshgrid(data.longitude, data.latitude);
% $$$ weights = 1./cosd(LAT(:));
% $$$ weights = mohsst5_remove_land(weights);
% $$$ 
% $$$ %
% $$$ % Means (uniform)
% $$$ %
% $$$ 
% $$$ t = 800;
% $$$ % t = 1200;
% $$$ 
% $$$ cl = max(max(abs(Y(:,t))), ...
% $$$          max(abs(GP1.F(:,t))));
% $$$ % $$$ cl = max(max(abs(Y(:,t))), ...
% $$$ % $$$          max(max(abs(GP1.F(:,t))), ...
% $$$ % $$$              max(abs(PCA1.F(:,t)))));
% $$$ 
% $$$ dir = '/home/jluttine/papers/gpkron/figures_aistats2012/';
% $$$ 
% $$$ fig();
% $$$ mohsst5_mapplot(GP1.F(:,t));
% $$$ colormap_redblue();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [-cl cl]);
% $$$ print('-depsc2', [dir 'fig_sst_gp_mean']);
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ mohsst5_mapplot(PCA1.F(:,t))
% $$$ % $$$ colormap_redblue();
% $$$ % $$$ polish(gca);
% $$$ % $$$ set(gca, 'clim', [-cl cl]);
% $$$ % $$$ print('-depsc2', [dir 'fig_sst_pca_mean']);
% $$$ 
% $$$ %
% $$$ % DATA
% $$$ %
% $$$ 
% $$$ fig()
% $$$ mohsst5_mapplot(Y(:,t))
% $$$ colormap_redblue();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [-cl cl]);
% $$$ print('-depsc2', [dir 'fig_sst_data']);
% $$$ 
% $$$ %
% $$$ % STD (uniform)
% $$$ %
% $$$ 
% $$$ burnin = 3;
% $$$ GPstd = sqrt(GP1.stdF(:,t).^2 + weights*mean(GP1.results.theta(end,burnin:end).^2));
% $$$ % $$$ PCAstd = sqrt(PCA1.stdF(:,t).^2 + 1./PCA1.results.Tau(:,t));
% $$$ 
% $$$ % $$$ fig()
% $$$ % $$$ plot(GPstd(:))
% $$$ % $$$ return
% $$$ 
% $$$ % $$$ fig()
% $$$ % $$$ %plot(GPstd(:))
% $$$ % $$$ plot(GP1.stdF(:))
% $$$ % $$$ min(GP1.stdF(:))
% $$$ % $$$ size(GP1.results.theta)
% $$$ % $$$ min(weights(:)*mean(GP1.results.theta(end,1000:end).^2))
% $$$ 
% $$$ 
% $$$ cl = max(GPstd(:));
% $$$ % $$$ cl = max( max(GPstd), max(PCAstd) );
% $$$ 
% $$$ fig()
% $$$ mohsst5_mapplot(GPstd);
% $$$ colormap_scale();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [0 cl]);
% $$$ print('-depsc2', [dir 'fig_sst_gp_std']);
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ mohsst5_mapplot(PCAstd);
% $$$ % $$$ colormap_scale();
% $$$ % $$$ polish(gca);
% $$$ % $$$ set(gca, 'clim', [0 cl]);
% $$$ % $$$ print('-depsc2', [dir 'fig_sst_pca_std']);
% $$$ 
% $$$ disp(datestr(data.time(t)));
% $$$ 
% $$$ % return
% $$$ 
% $$$ %return
% $$$ 
% $$$ %
% $$$ % Means (pattern)
% $$$ %
% $$$ 
% $$$ t = 1429;
% $$$ 
% $$$ cl = max(abs(GP2.F(:,t)));
% $$$ % $$$ cl = max( max(abs(GP2.F(:,t))), max(abs(PCA2.F(:,t))) );
% $$$ 
% $$$ % GP
% $$$ 
% $$$ fig()
% $$$ mohsst5_mapplot(GP2.F(:,t));
% $$$ colormap_redblue();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [-cl cl]);
% $$$ print('-depsc2', [dir 'fig_sst_pattern_gp_mean']);
% $$$ 
% $$$ % $$$ % PCA
% $$$ % $$$ 
% $$$ % $$$ figure
% $$$ % $$$ mohsst5_mapplot(PCA2.F(:,t))
% $$$ % $$$ colormap_redblue();
% $$$ % $$$ polish(gca);
% $$$ % $$$ set(gca, 'clim', [-cl cl]);
% $$$ % $$$ print('-depsc2', [dir 'fig_sst_pattern_pca_mean']);
% $$$ 
% $$$ % Data
% $$$ 
% $$$ fig()
% $$$ mohsst5_mapplot(Y(:,t))
% $$$ colormap_redblue();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [-cl cl]);
% $$$ print('-depsc2', [dir 'fig_sst_pattern_data']);
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ mohsst5_mapplot(PCA2.Ytrain(:,t))
% $$$ % $$$ colormap_redblue();
% $$$ % $$$ polish(gca);
% $$$ % $$$ set(gca, 'clim', [-cl cl]);
% $$$ % $$$ print('-depsc2', [dir 'fig_sst_pattern_train']);
% $$$ 
% $$$ disp(datestr(data.time(t)));
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ EGP = abs(Y-GP2.F);
% $$$ % $$$ EPCA = abs(Y-PCA2.F);
% $$$ % $$$ mohsst5_mapplot(EGP(:,t)-EPCA(:,t))
% $$$ % $$$ colormap_redblue();
% $$$ % $$$ polish(gca);
% $$$ % $$$ %set(gca, 'clim', [-cl cl]);
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ mohsst5_mapplot(E2(:,t))
% $$$ % $$$ colormap_scale();
% $$$ % $$$ polish(gca);
% $$$ %set(gca, 'clim', [-cl cl]);
% $$$ 
% $$$ figure()
% $$$ burnin = ceil(size(GP2.results.theta,2)/2);
% $$$ plot_scatterhist(GP2.results.theta(:,burnin:end)');
% $$$ 
% $$$ figure()
% $$$ subplot(1,2,1)
% $$$ semilogy(GP2.results.theta');
% $$$ subplot(1,2,2)
% $$$ plot(acorr(GP2.results.theta(:,burnin:end)'));
% $$$ 
% $$$ return
% $$$ 
% $$$ %
% $$$ % STD (pattern)
% $$$ %
% $$$ 
% $$$ GPstd = sqrt(GP2.stdF(:,t).^2 + weights*mean(GP2.results.theta(4,1000:end).^2));
% $$$ PCAstd = sqrt(PCA2.stdF(:,t).^2 + 1./PCA2.results.Tau(:,t));
% $$$ 
% $$$ cl = max( max(GPstd), max(PCAstd) );
% $$$ %cl = max( max(GP2.stdF(:,t)), max(PCA2.stdF(:,t)) );
% $$$ 
% $$$ figure
% $$$ mohsst5_mapplot(GPstd);
% $$$ colormap_scale();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [0 cl]);
% $$$ print('-depsc2', [dir 'fig_sst_pattern_gp_std']);
% $$$ 
% $$$ figure
% $$$ mohsst5_mapplot(PCAstd);
% $$$ colormap_scale();
% $$$ polish(gca);
% $$$ set(gca, 'clim', [0 cl]);
% $$$ print('-depsc2', [dir 'fig_sst_pattern_pca_std']);
% $$$ 
% $$$ return

function show_results(model, testset, filename)

dir = '/home/jluttine/papers/gpkron/figures/';
printfile = [dir 'fig_sst'];
if testset == 1
  teststring = 'uniform';
  t = 800;
else
  teststring = 'pattern';
  t = 1429;
end

if nargin >= 3
  % Load results
  dir = '/home/jluttine/matlab/publications/aistats2012/';
  Q = load([dir filename], 'F');
  % RMSE
  [rmse_train, rmse_test] = aistats2012_mohsst5_rmse(Q.F, testset);
  fprintf('train=%.4f and test=%.4f\n', rmse_train, rmse_test);
  F = Q.F(:,t);
  plot_map(F);
  printfile = sprintf('%s_%s_%s', printfile, teststring, model);
  print('-depsc2', printfile);
  title([teststring ' ' model]);
else
  [Ytrain, Ytest, Itrain, Itest] = aistats2012_mohsst5_sets(testset);
  F = Ytrain;
  F(Itest) = Ytest(Itest);
  plot_map(F(:,t));
  print('-depsc2', sprintf('%s_%s_data', printfile, teststring));
  title([teststring ' data']);
  plot_map(Ytrain(:,t));
  print('-depsc2', sprintf('%s_%s_train', printfile, teststring));
  title([teststring ' train']);
end


function plot_map(F)
% Plot reconstruction
cl = 3.3;
f = fig();
mohsst5_mapplot(F);
colormap_redblue();
polish(gca);
set(gca, 'clim', [-cl cl]);
% NaNs to grey
whitebg(f, [0.7 0.7 0.7]);

function polish(hax)
hcb = colorbar('SouthOutside');
set(hcb, 'Units', 'normalized');
set_figure_size(10,8);
set_colorbar_position(hcb, hax, 0.03, 0.07);

function f = fig()
f = figure('Position', [1 300 400 300]);
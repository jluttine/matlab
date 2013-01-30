function [norot, rot] = nc2010_plot_rppca_toyexperiment(n,m,d,ncomps)
%
% nc2010_plot_toyexperiment(model,n,m,d,ncomps,flatW)
%

model = 'rppca';

filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                    'nc2010_%s_toyexperiment_n=%d_m=%d_d=%d_ncomps=%d'], ...
                   model,n,m,d,ncomps); 

norot = load([filename,'_rotate=0.mat'], 'duration', 'cost', 'rmse', 'rmse_test');
rot = load([filename,'_rotate=1.mat'], 'duration', 'cost', 'rmse', 'rmse_test');
plot_results(norot, rot, n,m,d,ncomps);

if nargout < 1
  clear norot;
  clear rot;
end
return
% $$$ 
% $$$ norot = load('rotationpaper_experiment_n=100_m=10_d=3_ncomps=9_rotate=0.mat');
% $$$ rot = load('rotationpaper_experiment_n=100_m=10_d=3_ncomps=9_rotate=1.mat');
% $$$ plot_results(norot, rot, 100,10,3,9);
% $$$ 
% $$$ norot = load('rotationpaper_experiment_n=1000_m=50_d=10_ncomps=49_rotate=0.mat');
% $$$ rot = load('rotationpaper_experiment_n=1000_m=50_d=10_ncomps=49_rotate=1.mat');
% $$$ plot_results(norot, rot, 1000,50,10,49);
% $$$ 
% $$$ norot = load('rotationpaper_experiment_n=1000_m=100_d=10_ncomps=30_rotate=0.mat');
% $$$ rot = load('rotationpaper_experiment_n=1000_m=100_d=10_ncomps=30_rotate=1.mat');
% $$$ plot_results(norot, rot, 1000,100,10,30);
% $$$ 
% $$$ return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_results(norot, rot, n,m,d,ncomps)

% $$$ norot.cost(norot.cost==-inf) = -1e10;
% $$$ rot.cost(rot.cost==-inf) = -1e10;
suffix = sprintf(' (n=%d m=%d d=%d ncomps=%d)',n,m,d,ncomps);

figure
semilogx([norot.duration(:), rot.duration(:)], -[norot.cost(:), rot.cost(:)])
%title(['Negative loglikelihood', suffix]);
filename = sprintf('fig_cost_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 8,8]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),8,8])

figure
semilogx([norot.duration(:), rot.duration(:)], [norot.rmse(:), rot.rmse(:)])
%title(['Training set RMSE', suffix]);
filename = sprintf('fig_rmse_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 8,8]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),8,8])


figure
semilogx([norot.duration(:), rot.duration(:)], [norot.rmse_test(:), rot.rmse_test(:)])
%title(['Test set RMSE', suffix]);
filename = sprintf('fig_rmsetest_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 8,8]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),8,8])


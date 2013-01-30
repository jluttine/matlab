function figh = nc2010_plot_mlexperiment(flatW,ncomps,varargin)
%
%

% $$$ norot = load('rotationpaper_experiment_n=10_m=5_d=2_ncomps=4_rotate=0.mat');
% $$$ rot = load('rotationpaper_experiment_n=10_m=5_d=2_ncomps=4_rotate=1.mat');
% $$$ plot_results(norot, rot);
model = 'vbpca';
if flatW
  stringW = 'flatW';
  filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                      'nc2010_%s_mnistexperiment_%s_ncomps=%d'], ...
                     model,stringW,ncomps); 
else
  stringW = 'hierW';
  filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                      'nc2010_%s_mnistexperiment_%s_ncomps=%d'], ...
                     model,stringW,ncomps); 
end


norot = load([filename,'_rotate=0.mat'], 'duration', 'cost', 'rmse', 'rmse_test');
rot = load([filename,'_rotate=1.mat'], 'duration', 'cost', 'rmse', 'rmse_test');
figh = plot_results(norot, rot,ncomps, varargin{:});

% $$$ if nargout < 1
% $$$   clear norot;
% $$$   clear rot;
% $$$ end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figh = plot_results(norot, rot,ncomps, varargin)

opts = struct('xlim', [], ...
              'costlim', [], ...
              'rmselim', [], ...
              'rmsetestlim', [], ...
              'xtick', []);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end


% $$$ norot.cost(norot.cost==-inf) = -1e10;
% $$$ rot.cost(rot.cost==-inf) = -1e10;
%suffix = sprintf(' (n=%d m=%d d=%d ncomps=%d)',n,m,d,ncomps);

figure
figh = gcf;
semilogx(norot.duration(:), -norot.cost(:), 'k-')
hold on
semilogx(rot.duration(:), -rot.cost(:), 'k:', 'linewidth', 2)
%semilogx([norot.duration(:), rot.duration(:)], -[norot.cost(:), rot.cost(:)])
%title(['Negative loglikelihood', suffix]);
%filename = sprintf('fig_cost_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 8,8]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),8,8])
xlabel('time (s)')
if ~isempty(opts.xlim)
  set(gca, 'xlim', opts.xlim);
end
if ~isempty(opts.costlim)
  set(gca, 'ylim', opts.costlim);
end
if ~isempty(opts.xtick)
  set(gca, 'xtick', opts.xtick);
end

figure
figh = [figh gcf];
semilogx(norot.duration(:), norot.rmse(:), 'k-')
hold on
semilogx(rot.duration(:), rot.rmse(:), 'k:', 'linewidth', 2)
%semilogx([norot.duration(:), rot.duration(:)], [norot.rmse(:), rot.rmse(:)])
%title(['Training set RMSE', suffix]);
%filename = sprintf('fig_rmse_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 8,8]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),8,8])
xlabel('time (s)')
if ~isempty(opts.xlim)
  set(gca, 'xlim', opts.xlim);
end
if ~isempty(opts.rmselim)
  set(gca, 'ylim', opts.rmselim);
end
if ~isempty(opts.xtick)
  set(gca, 'xtick', opts.xtick);
end


figure
figh = [figh gcf];
semilogx(norot.duration(:), norot.rmse_test(:), 'k-')
hold on
semilogx(rot.duration(:), rot.rmse_test(:), 'k:', 'linewidth', 2)
%semilogx([norot.duration(:), rot.duration(:)], [norot.rmse_test(:), rot.rmse_test(:)])
%title(['Test set RMSE', suffix]);
%filename = sprintf('fig_rmsetest_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1:2), 8,8]);
pos = get(gcf, 'paperposition');
set(gcf, 'paperposition', [pos(1:2),8,8])
xlabel('time (s)')
if ~isempty(opts.xlim)
  set(gca, 'xlim', opts.xlim);
end
if ~isempty(opts.rmsetestlim)
  set(gca, 'ylim', opts.rmsetestlim);
end
if ~isempty(opts.xtick)
  set(gca, 'xtick', opts.xtick);
end

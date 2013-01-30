function figh = nc2010_plot_varianceexperiment(flatW,datatype,varargin)
%
% nc2010_plot_toyexperiment(model,n,m,d,ncomps,flatW)
%

debug = true

if debug
  debugstring = '';
else
  debugstring = 'debug';
end

n = 200;
m = 50;
d = 10;
ncomps = 30;


model = 'vbpca';

if flatW
  stringW = 'flatW';
else
  stringW = 'hierW';
end
filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                    'nc2010_%s_variance%sexperiment_%s_subspace=%s_n=%d_m=' ...
                    '%d_d=%d_ncomps=%d'], model,debugstring,stringW, ...
                   datatype,n,m,d,ncomps);
    

norot = load([filename,'_rotate=0.mat'], 'duration', 'cost', 'rmse', ...
             'rmse_test', 'w');
rot = load([filename,'_rotate=1.mat'], 'duration', 'cost', 'rmse', ...
           'rmse_test', 'w');
figh = plot_results(norot, rot, n,m,d,ncomps, varargin{:});

%norot.w
%rot.w

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figh = plot_results(norot, rot, n,m,d,ncomps, varargin)

opts = struct('xlim', [], ...
              'costlim', [], ...
              'rmselim', [], ...
              'rmsetestlim', []);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

% $$$ norot.cost(norot.cost==-inf) = -1e10;
% $$$ rot.cost(rot.cost==-inf) = -1e10;
suffix = sprintf(' (n=%d m=%d d=%d ncomps=%d)',n,m,d,ncomps);

if length(norot.duration) ~= length(rot.duration)
  ln = min(length(norot.duration), length(rot.duration));
  norot.duration = norot.duration(1:ln);
  norot.cost = norot.cost(1:ln);
  norot.rmse = norot.rmse(1:ln);
  norot.rmse_test = norot.rmse_test(1:ln);
  rot.cost = rot.cost(1:ln);
  rot.rmse = rot.rmse(1:ln);
  rot.rmse_test = rot.rmse_test(1:ln);
end

figure
figh = gcf;
norot.cost(norot.cost<-1e10) = nan;
rot.cost(rot.cost<-1e10) = nan;
semilogx(norot.duration(:), -norot.cost(:), 'k-')
hold on
semilogx(rot.duration(:), -rot.cost(:), 'k:', 'linewidth', 2)
%semilogx([norot.duration(:), rot.duration(:)], -[norot.cost(:), rot.cost(:)])
%title(['Negative loglikelihood', suffix]);
filename = sprintf('fig_cost_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
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

figure
figh = [figh gcf];
semilogx(norot.duration(:), norot.rmse(:), 'k-')
hold on
semilogx(rot.duration(:), rot.rmse(:), 'k:', 'linewidth', 2)
%semilogx([norot.duration(:), rot.duration(:)], [norot.rmse(:), rot.rmse(:)])
%title(['Training set RMSE', suffix]);
filename = sprintf('fig_rmse_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
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


figure
figh = [figh gcf];
semilogx(norot.duration(:), norot.rmse_test(:), 'k-')
hold on
semilogx(rot.duration(:), rot.rmse_test(:), 'k:', 'linewidth', 2)
%semilogx([norot.duration(:), rot.duration(:)], [norot.rmse_test(:), rot.rmse_test(:)])
%title(['Test set RMSE', suffix]);
filename = sprintf('fig_rmsetest_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
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



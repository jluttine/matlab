function nc2010_plot_compexperiment(flatW, datatype, convergence)

n = 200;
m = 50;
d = 10;

model = 'vbpca';
if flatW
  stringW = 'flatW';
else
  stringW = 'hierW';
end

% SELECT THE SEED
seed = 4132
datasets = 10
% Initialise variables
norottime = zeros(m,datasets);
norotloglike = zeros(m,datasets);
rottime = zeros(m,datasets);
rotloglike = zeros(m,datasets);
for ncomps = 1:m
  for ind = 1:datasets
    rotstring = 'rotate=0';
    norotfilename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                        'nc2010_%s_compexperiment_%s_subspace=%s_n=%d_m=' ...
                        '%d_d=%d_ncomps=%d_%s_seed=%d_datasetindex=%d.mat'], ...
                            model,stringW,datatype,n,m,d, ncomps, rotstring, ...
                            seed, ind);

    rotstring = 'rotate=1';
    rotfilename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                        'nc2010_%s_compexperiment_%s_subspace=%s_n=%d_m=' ...
                        '%d_d=%d_ncomps=%d_%s_seed=%d_datasetindex=%d.mat'], ...
                            model,stringW,datatype,n,m,d, ncomps, rotstring, ...
                            seed, ind);
    norot = load(norotfilename, 'duration', 'cost', 'rmse', 'rmse_test');
    rot = load(rotfilename, 'duration', 'cost', 'rmse', 'rmse_test');
    %plot_results(norot, rot,ncomps);
    
    % Find the point of convergence
    switch 2
     case 1 % by relative change in loglikelihood
      len = length(norot.cost);
      dcost = diff(norot.cost);
      dstep = abs(dcost ./ norot.cost(1:(len-1)));
      norotind = find(dstep < convergence, 1, 'first');
      len = length(rot.cost);
      dcost = diff(rot.cost);
      dstep = abs(dcost ./ rot.cost(1:(len-1)));
      rotind = find(dstep < convergence, 1, 'first');
     case 2 % being inside a threshold compared to converged loglikelihood
      convind = find(~isnan(norot.cost), 1, 'last');
      dif = abs((norot.cost - norot.cost(convind))/norot.cost(convind));
      norotind = find(dif < convergence, 1, 'first');
      convind = find(~isnan(rot.cost), 1, 'last');
      dif = abs((rot.cost - rot.cost(convind))/rot.cost(convind));
      rotind = find(dif < convergence, 1, 'first');
    end
    
    
    %  norotind = find(~isnan(norot.cost), 1, 'last');
    %  rotind = find(~isnan(rot.cost), 1, 'last');
    norottime(ncomps,ind) = norot.duration(norotind);
    norotloglike(ncomps,ind) = norot.cost(norotind);
    rottime(ncomps,ind) = rot.duration(rotind);
    rotloglike(ncomps,ind) = rot.cost(rotind);
  end
end

figure
norotmedian = prctile(norottime,50,2);
norot90 = prctile(norottime,100,2);
norot10 = prctile(norottime,0,2);
%norotmean = mean(norottime,2);
%rotmean = mean(rottime,2);
rotmedian = prctile(rottime,50,2);
rot90 = prctile(rottime,100,2);
rot10 = prctile(rottime,0,2);
semilogy(1:ncomps, norotmedian(:), 'r-');
hold on
x = [1:ncomps, ncomps:-1:1]';
y = [rot10(:); rot90(end:-1:1)];
c = [0.75 0.75 0.75];
fill(x,y,c,'EdgeColor',c);
%patch(x,y, [0.6 0.6 0.6]); 
% $$$ semilogy(1:ncomps, rot90(:), 'k:', 'linewidth', 1);
% $$$ semilogy(1:ncomps, rot10(:), 'k:', 'linewidth', 1);
semilogy(1:ncomps, norotmedian(:), 'r-');
semilogy(1:ncomps, norot90(:), 'r-');
semilogy(1:ncomps, norot10(:), 'r-');
semilogy(1:ncomps, rotmedian(:), 'k:', 'linewidth', 2);
%semilogy(1:ncomps, rot90(:), 'k:', 'linewidth', 1);
%semilogy(1:ncomps, rot10(:), 'k:', 'linewidth', 1);
%[norotloglike(:) - rotloglike(:)];
xlabel('number of components, D')
ylabel('time (seconds)')
xlim([1 m]);

return 

% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function plot_results(norot, rot,ncomps)
% $$$ 
% $$$ % $$$ norot.cost(norot.cost==-inf) = -1e10;
% $$$ % $$$ rot.cost(rot.cost==-inf) = -1e10;
% $$$ %suffix = sprintf(' (n=%d m=%d d=%d ncomps=%d)',n,m,d,ncomps);
% $$$ 
% $$$ figure
% $$$ subplot(1,3,1)
% $$$ semilogx([norot.duration(:), rot.duration(:)], -[norot.cost(:), rot.cost(:)])
% $$$ %title(['Negative loglikelihood', suffix]);
% $$$ %filename = sprintf('fig_cost_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
% $$$ % $$$ set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
% $$$ % $$$ pos = get(gcf, 'position');
% $$$ % $$$ set(gcf, 'position', [pos(1:2), 8,8]);
% $$$ % $$$ pos = get(gcf, 'paperposition');
% $$$ % $$$ set(gcf, 'paperposition', [pos(1:2),8,8])
% $$$ 
% $$$ %figure
% $$$ subplot(1,3,2)
% $$$ semilogx([norot.duration(:), rot.duration(:)], [norot.rmse(:), rot.rmse(:)])
% $$$ % $$$ %title(['Training set RMSE', suffix]);
% $$$ % $$$ %filename = sprintf('fig_rmse_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
% $$$ % $$$ set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
% $$$ % $$$ pos = get(gcf, 'position');
% $$$ % $$$ set(gcf, 'position', [pos(1:2), 8,8]);
% $$$ % $$$ pos = get(gcf, 'paperposition');
% $$$ % $$$ set(gcf, 'paperposition', [pos(1:2),8,8])
% $$$ 
% $$$ 
% $$$ %figure
% $$$ subplot(1,3,3)
% $$$ semilogx([norot.duration(:), rot.duration(:)], [norot.rmse_test(:), rot.rmse_test(:)])
% $$$ % $$$ %title(['Test set RMSE', suffix]);
% $$$ % $$$ %filename = sprintf('fig_rmsetest_n=%d_m=%d_d=%d_ncomps=%d', n,m,d,ncomps);
% $$$ % $$$ set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
% $$$ % $$$ pos = get(gcf, 'position');
% $$$ % $$$ set(gcf, 'position', [pos(1:2), 8,8]);
% $$$ % $$$ pos = get(gcf, 'paperposition');
% $$$ % $$$ set(gcf, 'paperposition', [pos(1:2),8,8])
% $$$ 
% $$$ 
function nc2010_publish_realexperiment

costlim = [3.6e5 4.1e5];
rmselim = [15 20];
rmsetestlim = [35.5 37.5];
xlim = [1e1 1e5];
xtick = 10.^(1:5);

filename = ['/home/jluttine/papers/neurocomputing2010/' ...
            'fig_realexperiment_mnist'];

figh = nc2010_plot_mnistexperiment(false, 50, 'xlim', xlim, 'costlim', ...
                                   costlim, 'rmselim', rmselim, 'rmsetestlim', ...
                                   rmsetestlim, 'xtick', xtick);

for ind = 1:3
  set(figh(ind), 'units', 'centimeters', 'paperunits', 'centimeters');
  pos = get(figh(ind), 'position');
  set(figh(ind), 'position', [pos(1:2), 8,6]);
  pos = get(figh(ind), 'paperposition');
  set(figh(ind), 'paperposition', [pos(1:2),8,6])
  
  switch ind
   case 1
    print(figh(ind), '-depsc2', [filename, '_cost']);
   case 2
    print(figh(ind), '-depsc2', [filename, '_rmse']);
   case 3
    print(figh(ind), '-depsc2', [filename, '_rmsetest']);
  end
end

costlim = [1e5 5e5];
rmselim = [0.5 1];
rmsetestlim = [0.88 0.96];
xlim = [1e2 1e4];
xtick = [];

filename = ['/home/jluttine/papers/neurocomputing2010/' ...
            'fig_realexperiment_movielens'];

figh = nc2010_plot_mlexperiment(false, 100, 'xlim', xlim, 'costlim', costlim, ...
                                'rmselim', rmselim, 'rmsetestlim', ...
                                rmsetestlim, 'xtick', xtick);
for ind = 1:3
  set(figh(ind), 'units', 'centimeters', 'paperunits', 'centimeters');
  pos = get(figh(ind), 'position');
  set(figh(ind), 'position', [pos(1:2), 8,6]);
  pos = get(figh(ind), 'paperposition');
  set(figh(ind), 'paperposition', [pos(1:2),8,6])
  
  switch ind
   case 1
    print(figh(ind), '-depsc2', [filename, '_cost']);
   case 2
    print(figh(ind), '-depsc2', [filename, '_rmse']);
   case 3
    print(figh(ind), '-depsc2', [filename, '_rmsetest']);
  end
end


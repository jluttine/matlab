function nc2010_publish_toyexperiment

datatype = 'weak';

for flatW = [true false]
  
  if flatW 
    stringW = 'flatW';
    costlim = [];
    rmselim = [];
    rmsetestlim = [1.1 1.7];
    xlim = [1e0 1e3];
  else
    stringW = 'hierW';
    costlim = [];
    rmselim = [];
    rmsetestlim = [1.1 1.7];
    xlim = [1e0 1e3];
  end
  
  filename = sprintf(['/home/jluttine/papers/neurocomputing2010/' ...
                      'fig_toyexperiment_%s_subspace=%s'], stringW, datatype);

  figh = nc2010_plot_toyexperiment(200,50,10,30, flatW, datatype, 'xlim', ...
                                   xlim, 'costlim', costlim, 'rmselim', ...
                                   rmselim, 'rmsetestlim', rmsetestlim);
  
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

end

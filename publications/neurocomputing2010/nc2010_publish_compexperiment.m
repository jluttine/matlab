function nc2010_publish_compexperiment

threshold = 1e-3


for flatW = [true false]
  for datatype = {'no', 'weak', 'strong'}
    
    if flatW 
      stringW = 'flatW';
    else
      stringW = 'hierW';
    end
      
    filename = sprintf(['/home/jluttine/papers/neurocomputing2010/' ...
                        'fig_compexperiment_severalruns_%s_subspace=%s'], ...
                       stringW, datatype{:})

    nc2010_plot_compexperiment(flatW, datatype{:}, threshold);
    
    set(gcf, 'units', 'centimeters', 'paperunits', 'centimeters');
    pos = get(gcf, 'position');
    set(gcf, 'position', [pos(1:2), 8,6]);
    pos = get(gcf, 'paperposition');
    set(gcf, 'paperposition', [pos(1:2),8,6])

    print(gcf, '-depsc2', filename);
  end
end

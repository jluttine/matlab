
function metoffice_experiment_comparison(remove_val_set, anomalies_set)

% Runs GPFA and GPloc with different initializations and parameter
% values for a simulated data.  Then, validation errors and test errors
% can be compared.

seed = floor(10000*rand());

dataset = 'hadsst2d1';

maxiter = 500;
N_samples = 2000;

if nargin < 1
  remove_val_set = [2,4];
end
if nargin < 2
  anomalies_set = [1,0];
end

for remove_val = remove_val_set
  for anomalies = anomalies_set
    
    % Load data
    data = metoffice_get_data(dataset, anomalies, remove_val);
    % Form the data matrix
    Y = data.data;
    times = data.time;
    coordinates = data.coordinates;


    % Local GP
    if anomalies
      randn('state', seed);
      rand('state', seed);
      res = metoffice_experiment_gpkron(Y, ...
                                        times, ...
                                        coordinates, ...
                                        N_samples, ...
                                        '', ...
                                        '');
      Yh = res.F;
      clear res
      performance_measure(Yh, ...
                          dataset, ...
                          anomalies, ...
                          remove_val, ...
                          'gpkron', ...
                          ['D=1_seed=', num2str(seed)]);
    end
    
    % GPFA
    for D = [20, 40, 60, 80]
      if anomalies
        comps_temporal = [0 10 D-10];
        comps_spatial =  [0 D  0];
      else
        comps_temporal = [10 10   D-20];
        comps_spatial =  [10 D-10 0];
      end
      randn('state', seed);
      rand('state', seed);
      Q = metoffice_experiment_gpfa(Y, ...
                                    comps_temporal, ...
                                    comps_spatial, ...
                                    times, ...
                                    coordinates, ...
                                    maxiter, ...
                                    '');
      Yh = Q.W'*Q.X;
      clear Q
      performance_measure(Yh, ...
                          dataset, ...
                          anomalies, ...
                          remove_val, ...
                          'gpfa', ...
                          sprintf('D=%d_seed=%d', ...
                                  sum(comps_temporal), ...
                                  seed));
    end
    
  end
end

                   



function perfmeas = performance_measure(Yrec, dataset, anomalies, remove_val, ...
                                        method, suffix)

%global datapath
%datapath = '/share/climate/data/UK_Met_Office';
%addpath /home/alexilin/matlab/metoffice
Yrec = metoffice_add_land( Yrec, dataset );
Yrec = metoffice_add_climatology( Yrec, dataset, anomalies );
% $$$ Yrec = add_land( Yrec, dataset );
% $$$ Yrec = add_climatology( Yrec, dataset, anomalies );
perfmeas = metoffice_compute_rmse( Yrec, dataset, remove_val );
perfmeas.setup = struct('dataset', dataset, ...
                        'anomalies', anomalies, ...
                        'validation', remove_val, ...
                        'method', method, ...
                        'other', suffix);

folder = '/home/jluttine/matlab/datasets/metoffice/comparison_results';
filename = sprintf('%s/metoffice_comparison_dataset=%s_validation=%d_method=%s_anomalies=%d_%s', ...
                   folder, ...
                   dataset, ...
                   remove_val, ...
                   method, ...
                   anomalies, ...
                   suffix);
save(filename, 'perfmeas');

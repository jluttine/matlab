% Run experiments for simulated datasets.
% Usage: Q = metoffice_experiment_pca('hadsst2d1', 1, 40, 4)

function [ Q, perfmeas ] = metoffice_experiment_pca(dataset, anomalies, D, remove_val)

%
% Load data
%

if nargin < 4 || isempty(remove_val)
    remove_val = 1;
end
[data,folder,maskfile] = metoffice_getdata(dataset, anomalies, remove_val);

% Form the data matrix
Y = data.data;
[M,N] = size(Y);
Obs = ~isnan(Y);

%
% PCA inference
%

% Number of components
if nargin < 3 || isempty(D)
    D = 120;
end

% Filename for saving the results
%folder = [folder '/pca'];
mkdir(folder);
filename = sprintf('%s/%s/pca_D=%d_anomalies=%d_remval=%02d_%s', ...
                   folder, ...
                   dataset, ...
                   D, ...
                   anomalies, ...
                   remove_val, ...
                   datestr(now,'yyyymmdd'))

% PCA module for X (one constant component for modeling bias)
prior.mu = [1; zeros(D-1,1)];
prior.CovX = diag([1e-6; ones(D-1,1)]);
X_module = factor_module_iid('prior', prior);

% ARD module for W
W_module = factor_module_ard();

% Isotropic noise (precisions weighted proportionally to grid size)
[LON,LAT] = meshgrid(data.lon, data.lat);
weights = cosd(LAT(:));
weights = metoffice_remove_bins(weights,maskfile);
weights = repmat(weights, [1, N]);
noise_module = noise_module_isotropic('init', struct('tau', 10), ...
                                      'weights', weights);

% Run VB PCA
Q = vbfa(D, Y, W_module, X_module, noise_module, ...
         'maxiter', 2, ...
         'rotate', true, ...
         'autosavefile', filename, ...
         'autosave', [1 20:20:2000]);

if nargout == 2
    % Reconstruct
    Yrec = Q.W'*Q.X;
    
    % Remove the autosave file
    delete( [ filename '.mat' ] )

    addpath /home/alexilin/matlab/metoffice
    Yrec = add_land( Yrec, dataset );
    Yrec = add_climatology( Yrec, dataset, anomalies );
    perfmeas = compute_rmse( Yrec, dataset, remove_val );

else
    % Reconstruct
    Yrec = Q.W'*Q.X;
    
    % Some performance measures
    fprintf('Weighted training RMSE of the reconstruction: %f\n',  ...
            rmsew(Y(Obs)-Yrec(Obs),weights(Obs)));

end

% Save the results
%save(filename, '-struct', 'Q');

function Q = testbed_experiment_vbrfa(d, common_nu, common_tau, maxiter, varargin)

randn('state', 666);
rand('state', 666);

% Default parameters
options = struct( ...
    'autosavetime', 100, ...
    'rotate',       1, ...
    'update_nu',    [2, 3, 4, 5:5:maxiter], ...
    'update_alpha',     1, ...    % don't update before starting rotations
    'robustness',   'independent-t', ...
    'prior',        []);

% Parse arguments
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Get data
disp('Loading data..');
data = testbed_loaddata();
data = testbed_preprocess(data);

filename = sprintf(['/share/climate/jluttine/testbed/' ...
                    'testbed_results_d=%d_robustness=%s_common-nu=' ...
                    '%d_common-tau=%d'], d, options.robustness, ...
                   common_nu, common_tau);

% Remove 48th station?
remove48 = false;
if remove48
  filename = sprintf('%s_remove48', filename);
  data = testbed_remove_stations(data, 48);
  disp('NOTE: Removed station 48.');
end

% Remove stations 20-21?
remove2021 = false;
if remove2021
  filename = sprintf('%s_remove2021', filename);
  data = testbed_remove_stations(data, 20:21);
  disp('NOTE: Removed stations 20 and 21.');
end

% Use prior?
use_prior = false;
if use_prior
  filename = sprintf('%s_useprior', filename);
  options.prior.mumu = 0;
  options.prior.taumu = 1/(4^2);
end

% Show the options
options.common_nu = common_nu;
options.common_tau = common_tau;
options.user_data = data;
options.init.nu = 1;
options.init.tau = 1e-0;
options.maxiter = maxiter;

% Initialize with 1D results?
init1D = false;
if init1D
  options2.prior = options.prior;
  options2.common_nu = options.common_nu;
  options2.common_tau = options.common_tau;
  options2.maxiter = 10;
  options2.update_nu = 2:options2.maxiter;%1:options.maxiter;
  options2.init.nu = 1;
  options2.init.tau = options.init.tau;
  options2.user_data = options.user_data;
  
  Q = vbrfa(data.observations, 1, options2);
  
  % Use the results to initialize 
  options.init.U = Q.U;
  options.init.nu = Q.nu;
  options.init.tau = Q.tau;
  filename = sprintf('%s_init1d', filename);
end

filename = sprintf('%s_date=%s', filename, datestr(now, 'yyyymmdd'));
options.autosavefile = filename;

fprintf('Parameters: d=%d, common_nu=%d, common_tau=%d, maxiter=%d\n', d, ...
        common_nu, common_tau, maxiter);

fprintf('Results will be saved in %s.mat\n', filename);

if nargin < 2
  maxiter = 1000;
end

% Show the options
options

disp('Running variational Bayesian robust factor analysis..');
Q = vbrfa(data.observations, d, options);

save(filename, '-struct', 'Q');



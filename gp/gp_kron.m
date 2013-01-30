% TODO:
%
% MAYBE THE HYPERPARAMETER SAMPLING COULD BE MADE MORE EFFICIENT BY USING
% SIMULATED ANNEALING?!?!

function res = gp_kron(Y, covfunc1, covfunc2, N_samples, theta_init, ...
                       logprior_theta, dlogprior_theta, samplefunc, varargin)

%
% DATA
%

% Parse options
options = struct('hyperparameters', 'ss', ...
                 'filename_samples', [], ...
                 'noise_scale1', [], ...
                 'noise_scale2', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

[N1,N2] = size(Y);

Imv = isnan(Y);

%
% INFERENCE
%

switch 'kron6'
% $$$  case 'kron1'
% $$$ 
% $$$   a = 1e-3;
% $$$   b = 1e-3;
% $$$   logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
% $$$   dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);
% $$$ 
% $$$   % Initial guess for covariance parameters
% $$$   theta_init = [0.3         ... % overall magnitude
% $$$                 0.5*theta1  ... % length scale
% $$$                 1.3*sqrt(s) ... % noise magnitude
% $$$                 0.5*theta2  ... % length scale
% $$$                 1.3*sqrt(s)]';  % noise magnitude
% $$$   
% $$$   samplefunc = get_sample_function2(numel(theta_init), N_samples, burnin);
% $$$   
% $$$   noise_weights = linspace(1,2,N1)';
% $$$ 
% $$$   [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
% $$$       gp_init_kron1(covfunc1, ...
% $$$                     solver_ldlchol(), ...
% $$$                     covfunc2, ...
% $$$                     solver_ldlchol(), ...
% $$$                     logprior_theta, ...
% $$$                     dlogprior_theta, ...
% $$$                     samplefunc, ...
% $$$                     'noise_scale1', noise_weights);
% $$$ 
% $$$  case 'kron2'
% $$$    
% $$$   a = 1e-3;
% $$$   b = 1e-3;
% $$$   logprior_theta = @(theta) sum(gamma_logpdf(theta, a, b));
% $$$   dlogprior_theta = @(theta) gamma_dlogpdf(theta, a, b);
% $$$ 
% $$$   % Initial guess for covariance parameters
% $$$   theta_init = [2.0        ... % total magnitude
% $$$                 2.0*theta1 ... % 1) length scale
% $$$                 0.5*theta2 ... % 2) length scale
% $$$                 2.0]';         % noise magnitude
% $$$   
% $$$   samplefunc = get_sample_function2(numel(theta_init), N_samples, burnin);
% $$$ 
% $$$ % $$$   [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
% $$$ % $$$       gp_init_kron2(covfunc1, ...
% $$$ % $$$                     solver_ldlchol(), ...
% $$$ % $$$                     covfunc2, ...
% $$$ % $$$                     solver_ldlchol(), ...
% $$$ % $$$                     logprior_theta, ...
% $$$ % $$$                     dlogprior_theta, ...
% $$$ % $$$                     'samplefunc', samplefunc, ...
% $$$ % $$$                     'rand_y', 'gibbs', ...
% $$$ % $$$                     'likelihood', 'conditioned');
% $$$ 
% $$$   [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
% $$$       gp_init_kron2(covfunc1, ...
% $$$                     solver_ldlchol(), ...
% $$$                     covfunc2, ...
% $$$                     solver_ldlchol(), ...
% $$$                     logprior_theta, ...
% $$$                     dlogprior_theta, ...
% $$$                     'samplefunc', samplefunc, ...
% $$$                     'rand_y', 'pcg', ...
% $$$                     'likelihood', 'whitened_approximate');

 case 'kron4'

  if isnumeric(samplefunc)
    burnin = samplefunc; %floor(N_samples/2);
    samplefunc = get_sample_function2(numel(theta_init), N_samples, burnin);
  end
  
  [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
      gp_init_kron4(covfunc1, ...
                    solver_ldlchol(), ...
                    covfunc2, ...
                    solver_ldlchol(), ...
                    logprior_theta, ...
                    dlogprior_theta, ...
                    'samplefunc', samplefunc, ...
                    'rand_y', 'pcg', ...
                    'noise_scale', options.noise_scale, ...
                    'likelihood', 'whitened_prior');
  
 case 'kron5'
  
  if isnumeric(samplefunc)
    burnin = samplefunc; %floor(N_samples/2);
    samplefunc = get_sample_function2(numel(theta_init), N_samples, burnin);
  end
  
  if iscell(covfunc1)
    covfuncs = {covfunc1{:}; covfunc2{:}};
  else
    covfuncs = {covfunc1; covfunc2};
  end
  solvers = cell(size(covfuncs));
  solvers(:) = {solver_ldlchol()};
  
  [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
      gp_init_kron5(covfuncs, ...
                    solvers, ...
                    logprior_theta, ...
                    dlogprior_theta, ...
                    'samplefunc', samplefunc, ...
                    'rand_y', 'pcg', ...
                    'noise_scale', options.noise_scale, ...
                    'likelihood', 'whitened_prior');

 case 'kron6'
  
  if isnumeric(samplefunc)
    burnin = samplefunc; %floor(N_samples/2);
    samplefunc = get_sample_function2(numel(theta_init), ...
                                      N_samples, ...
                                      burnin, ...
                                      'filename_samples', options.filename_samples);
  end
  
  if iscell(covfunc1)
    covfuncs = {covfunc1{:}; covfunc2{:}};
  else
    covfuncs = {covfunc1; covfunc2};
  end
  solvers = cell(size(covfuncs));
  solvers(:) = {solver_ldlchol()};
  
  [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
      gp_init_kron6(covfuncs, ...
                    solvers, ...
                    logprior_theta, ...
                    dlogprior_theta, ...
                    'samplefunc', samplefunc, ...
                    'rand_y', 'pcg', ...
                    'noise_scale1', options.noise_scale1, ...
                    'noise_scale2', options.noise_scale2, ...
                    'likelihood', 'whitened_prior');

end

% Transform to log-scale
func_theta = @(logtheta, varargin) func_theta_transformed(logtheta, ...
                                                  exp(logtheta), ...
                                                  func_theta, ...
                                                  varargin{:});
get_logpdf = @(f_theta) get_logpdf_transformed(f_theta, ...
                                               get_logpdf, ...
                                               sum(f_theta.theta_transformed));
get_dlogpdf = @(df_theta) get_dlogpdf_transformed(df_theta, ...
                                                  get_dlogpdf, ...
                                                  diag(exp(df_theta.theta_transformed)), ...
                                                  ones(size(df_theta.theta_transformed)));
theta_init = log(theta_init);

% Posterior sampling for covariance parameters
switch options.hyperparameters
  
% $$$  case 'debug'
% $$$ 
% $$$   [rand_theta, f_theta] = mcmc_init_reflective(theta_init, ...
% $$$                                                get_logpdf, ...
% $$$                                                get_dlogpdf, ...
% $$$                                                0.001, ...
% $$$                                                100, ...
% $$$                                                'type', 'outside', ...
% $$$                                                'fx', func_theta);
% $$$   
% $$$  case 'mh'
% $$$   % Metropolis-Hastings
% $$$ 
% $$$   % Gaussian proposals
% $$$   stepsize = 30e-3;
% $$$   rand_proposal = @(logtheta) normrnd(logtheta, stepsize);
% $$$   logpdf_proposal = @(logtheta1, logtheta2) ...
% $$$       sum(normal_logpdf(logtheta1,      ...
% $$$                         logtheta2, ...
% $$$                         stepsize));
% $$$ 
% $$$   [rand_theta, f_theta] = mcmc_init_metropolishastings(theta_init,   ...
% $$$                                                     get_logpdf,      ...
% $$$                                                     rand_proposal,   ...
% $$$                                                     logpdf_proposal, ...
% $$$                                                     func_theta);
% $$$ 
% $$$   
 case 'hmc'
  % Hamiltonian/Hybrid Monte Carlo
  
  [rand_theta, f_theta] = mcmc_init_hamiltonian(theta_init, ...
                                                get_logpdf, ...
                                                get_dlogpdf, ...
                                                1e-2, ...
                                                10, ...
                                                func_theta);
% $$$  case 'reflective'
% $$$ 
% $$$   [rand_theta, f_theta] = mcmc_init_reflective(theta_init, ...
% $$$                                                get_logpdf, ...
% $$$                                                get_dlogpdf, ...
% $$$                                                0.001, ...
% $$$                                                100, ...
% $$$                                                'type', 'inside', ...
% $$$                                                'fx', func_theta);
% $$$   
 case 'ss'
  % Slice sampling (do sampling in log-scale)
  %
  % NOTE: Slice sampling can be very bad if initialized at very low
  % density values. Sometimes it might be surprising that a "good"
  % initial guess happens to have very low density so the algorithm
  % won't work in practice.
  
  [rand_theta, f_theta] = mcmc_init_slicesampling(theta_init, ...
                                                  get_logpdf, ...
                                                  'fx', func_theta);
  
% $$$  case 'ml'
% $$$   % Maximum likelihood
% $$$   rand_theta = init_maximumlikelihood();
  
 case 'fixed'
  rand_theta = @fixed_theta
  f_theta = func_theta(theta_init);
  
 case default
  error('Unknown sampling scheme for theta');

end

rand_model = get_rand_model(rand_theta, f_theta);

% Gibbs sampling
t = cputime();
gibbs(Y,Imv, N_samples, rand_model);


%
% ANALYSIS
%

% Check the results
res = samplefunc();

    function [x, fx] = fixed_theta(varargin)
    x = theta_init;
    fx = f_theta;
    end

end

function S = solver_ldlchol()

S.decompose = @decompose;
S.linsolve = @linsolve;
S.logdet = @logdet;
S.inv = @inv;
S.squareroot = @squareroot;
  
  function s = decompose(K)
  [s.LD,s.p] = ldlchol(K);
  end
  
  function x = linsolve(s, y)
  x = linsolve_ldlchol(s.LD,y);
  end
  
  function ldet = logdet(s)
  ldet = logdet_ldlchol(s.LD);
  end
  
  function A = inv(s)
  A = spinv_ldlchol(s.LD);
  end
  
  function L = squareroot(s)
  L = ldlchol2lchol(s.LD);
  end

  %
  % WARNING: If you use permutation q as in [L,p,q]=ldlchol, note that
  % whitened loglikelihoods are messed up because the latent space
  % representation varies with the permutation.. Thus, you must solve the
  % Cholesky without using permutations for the squareroot.
  %

% $$$   function s = decompose(K)
% $$$   %s.q = 1:length(K);
% $$$   s.K = K;
% $$$   [s.LD,s.p,s.q] = ldlchol(K);
% $$$   end
% $$$   
% $$$   function x = linsolve(s, y)
% $$$   x(s.q,:) = linsolve_ldlchol(s.LD,y(s.q,:));
% $$$   end
% $$$   
% $$$   function ldet = logdet(s)
% $$$   ldet = logdet_ldlchol(s.LD);
% $$$   end
% $$$   
% $$$   function A = inv(s)
% $$$   error('Check this (q)')
% $$$   A(s.q,s.q) = spinv_ldlchol(s.LD);
% $$$   end
% $$$   
% $$$   function L = squareroot(s)
% $$$   L = lchol(s.K);
% $$$   end

end


function [f_theta, df_theta] = func_theta_transformed(theta_transformed, ...
                                                  theta, func_theta, varargin)
if nargout <= 1
  f_theta = func_theta(theta, varargin{:});
  f_theta.theta_transformed = theta_transformed;
else
  [f_theta, df_theta] = func_theta(theta, varargin{:});
  f_theta.theta_transformed = theta_transformed;
  df_theta.theta_transformed = theta_transformed;
end
end

function logpdf_y = get_logpdf_transformed(fy, get_logpdf, logjacobian)
logpdf = get_logpdf(fy);
logpdf_y = @logpdf_transformed;
  function lpdf = logpdf_transformed(varargin)
  lpdf = logpdf(varargin{:}) + logjacobian;
  end
end

function dlogpdf_y = get_dlogpdf_transformed(dfy, get_dlogpdf, Jacobian, ...
                                                  dlogjacobian)
dlogpdf = get_dlogpdf(dfy);
dlogpdf_y = @dlogpdf_transformed;
  function dlpdf = dlogpdf_transformed(varargin)
  dlpdf = dlogpdf(varargin{:});
  dlpdf = Jacobian*dlpdf + dlogjacobian;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samplefunc = get_sample_function2(D_theta, N, burnin, varargin)
% Parse options
options = struct('filename_samples', [], ...
                 'noise_scale2', []);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

results.Y = [];
results.F = 0;
results.FF = 0;
results.F2 = 0;
results.sumF = 0;
results.sumF2 = 0;
results.theta = zeros(D_theta, N);
samplefunc = @process_sample;
n = 1;
  function res = process_sample(Y, F, theta, F_mean)
  if nargin >= 1
    results.Y = Y;
    if  n > burnin
      results.F = (F + (n-burnin-1)*results.F) / (n-burnin);
      results.FF = (F.*F + (n-burnin-1)*results.FF) / (n-burnin);
      results.F2 = results.FF;
      if size(F,3) > 1
        sumF = sum(F,3);
        results.sumF = (sumF + (n-burnin-1)*results.sumF) / (n-burnin);
        results.sumF2 = (sumF.*sumF + (n-burnin-1)*results.sumF2) / (n-burnin);
      end
    end
    results.theta(:,n) = theta(:);
    if ~isempty(options.filename_samples)
      padded_zeros = ceil(log10(N+1));
      fmt = sprintf('_%%0%dd', padded_zeros);
      filename = [options.filename_samples, sprintf(fmt, n)];
      save(filename, 'F')
    end
  end
  if nargout >= 1
    res = results;
  end
  n = n + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gibbs(Y, Imv, N, rand_model, samplefunc)

Y(Imv) = 0;

for n=1:N
  
  t = cputime();
  Y = rand_model(Y,Imv);
  dt = cputime() - t;
  fprintf('Iteration step %d done. (%.3f seconds)\n', n, dt)
  
end

end





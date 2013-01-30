function Q = vbrfa2011_testbed_experiment(model, D, varargin)

options = struct('maxiter', 500);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

%
% Load the training data
%

%data_train = testbed_preprocess(testbed_loaddata());
data_train = load(['/share/climate/jluttine/testbed/' ...
                   'testbed_vbrfa2011_traindata']);

%
% Construct the model
%

% ARD for W and isotropic Gaussian noise
% $$$ W_module = ...
% $$$     factor_module_ard('noise_module', ...
% $$$                       noise_module_diagonal('init', struct('tau', 1), ...
% $$$                                             'separation', 'rows'));
W_module = factor_module_ard('noise_module', ...
                             noise_module_isotropic('init', ...
                                                  struct('tau', 0.01)));

% IID for X (do not model the bias term)
X_module = factor_module_iid();
% $$$ X_module = factor_module_iid('prior', struct('mu', [1; zeros(D,1)], ...
% $$$                                              'CovX', diag([1e-6; ones(D,1)])));



switch model

 case 1 % GAUSSIAN
  noise = noise_module_fixed(1);
  id = 'gaussian';

 case 2 % MULTIVARIATE T
  noise = noise_module_multivariate_t('update_nu', [5], ...
                                      'init', struct('nu', 0.1));
% $$$   noise = noise_module_product(...
% $$$       noise_module_multivariate_t('update_nu', [2:5, 10:5:1000], ...
% $$$                                   'init', struct('nu', 1)), ...
% $$$       isotropic_noise);
  id = 'multi-t';

 case 3 % INDEPENDENT T
  noise = noise_module_independent_t('update_nu', [5], ...
                                     'nu', 'separate_rows', ...
                                     'init', struct('nu', 0.1));
% $$$   noise = noise_module_product(...
% $$$       noise_module_independent_t('update_nu', [2:5 10:5:1000], ...
% $$$                                  'nu', 'separate_rows', ...
% $$$                                  'init', struct('nu', 1)), ...
% $$$       isotropic_noise);
  id = 'ind-t';

 case 4 % LAPLACE
  noise = noise_module_laplace();
% $$$   noise = noise_module_product(noise_module_laplace(), ...
% $$$                                isotropic_noise);
  id = 'laplace';
  
 otherwise
  error('Unknown model requested')
  
end

%
% Run the method
%

% Autosave-file
folder = '/share/climate/jluttine/testbed';
filename = sprintf('%s/testbed_results_%s_D=%d_%s', ...
                   folder, ...
                   id, ...
                   D, ...
                   datestr(now, 'yyyymmdd'));

% General function call for running the method
Q = vbfa(D, ...
         data_train.observations, ...
         W_module, ...
         X_module, ...
         noise, ...
         'maxiter', options.maxiter, ...
         'rotate', [1:9, 10:5:1000], ...
         'rotation_maxiter', 1000, ...
         'autosave', [1 5:5:1000], ...
         'autosavefile', filename); 

fprintf('Saving results to %s...', filename);
save(filename, '-struct', 'Q');
fprintf(' done.\n');

%
% Compute predictive RMSE
%

data_test = load(['/share/climate/jluttine/testbed/' ...
                  'testbed_vbrfa2011_testdata']);

Yh = Q.W'*Q.X;

Itrain = ~isnan(data_train.observations);
Itest = ~isnan(data_test.observations);
fprintf('RMSE-train: %.4f\n', ...
        rmse(data_train.observations(Itrain)-Yh(Itrain)));
fprintf('RMSE-test:  %.4f\n', ...
        rmse(data_test.observations(Itest)-Yh(Itest)));

%
% Compute predictive log-likelihood
%

disp('Computing predictive log-likelihood..')

vbrfa2011_testbed_measure(model, Q);

% $$$ 
% $$$ loglike = 0;
% $$$ samples = 100;
% $$$ [M,N] = size(Yh);
% $$$ Ytest = data_test.observations;
% $$$ for i=samples
% $$$   
% $$$   % Sample a noise level
% $$$   tau = gamrnd(Q.W_struct.rho_struct.a_tau, 1/Q.W_struct.rho_struct.b_tau);
% $$$ 
% $$$   % Sample a reconstruction
% $$$   for m=1:M
% $$$     %W(:,m) = mvnrnd(Q.W(:,m), 1/tau * Q.CovW(:,:,m));
% $$$     W(:,m) = gaussian_rand(Q.W(:,m), 1/tau * Q.CovW(:,:,m));
% $$$   end
% $$$   for n=1:N
% $$$ % $$$     X(:,n) = mvnrnd(Q.X(:,n), Q.CovX(:,:,n));
% $$$     X(:,n) = gaussian_rand(Q.X(:,n), Q.CovX(:,:,n));
% $$$   end
% $$$ % $$$   W = mvnrnd(Q.W', 1/tau * Q.CovW)';
% $$$ % $$$   X = mvnrnd(Q.X', Q.CovX)';
% $$$   F = W'*X;
% $$$ 
% $$$   % Log-likelihood
% $$$   switch model
% $$$    case 1 % GAUSSIAN
% $$$     Z = Ytest(Itest) - F(Itest);
% $$$     loglike = loglike + gaussian_logpdf(tau*(Z'*Z), ...
% $$$                                         0, ...
% $$$                                         0, ...
% $$$                                         -numel(Z)*log(tau), ...
% $$$                                         numel(Z));
% $$$    case 2 % MULTIVARIATE T
% $$$     alpha = Q.Tau_struct.a_u;
% $$$     beta = Q.Tau_struct.b_u;
% $$$     tau = tau * (alpha./beta);
% $$$     nu = 2*alpha;
% $$$     Z = Ytest - F;
% $$$     Z(~Itest) = 0;
% $$$     Z2 = tau .* dot(Z,Z,1);
% $$$     Nmv = sum(Itest,1)>0;
% $$$     loglike = loglike + sum(t_logpdf(Z2(Nmv), ...
% $$$                                      -log(tau(Nmv)), ...
% $$$                                      nu(Nmv), ...
% $$$                                      sum(Itest(:,Nmv),1)));
% $$$     
% $$$    case 3 % INDEPENDENT T
% $$$     [I,J] = find(Itest);
% $$$     Z2 = tau * (Ytest(Itest) - F(Itest)).^2;
% $$$     loglike = loglike + sum(t_logpdf(Z2, ...
% $$$                                      -log(tau), ...
% $$$                                      Q.Tau_struct.nu(I), ...
% $$$                                      1));
% $$$     
% $$$    case 4 % LAPLACE
% $$$     loglike = loglike + sum(laplace_logpdf(Ytest(Itest), ...
% $$$                                            F(Itest), ...
% $$$                                            sqrt(tau)));
% $$$   end
% $$$ 
% $$$ end
% $$$ loglike = loglike / samples;
% $$$ fprintf('Predictive log-likelihood: %.4e\n', loglike);

if nargout < 1
  clear Q;
end

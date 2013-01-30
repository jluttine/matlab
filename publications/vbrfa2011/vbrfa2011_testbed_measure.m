
function vbrfa2011_testbed_measure(model, Q)

if nargin < 2 || isempty(Q)
  % Load default results
  disp('Loading default results..')
  files = {'testbed_results_gaussian_D=30_20110617', ...
           'testbed_results_multi-t_D=30_20110617', ...
           'testbed_results_ind-t_D=30_20110617', ...
           'testbed_results_laplace_D=30_20110617'};
  Q = load(files{model}, 'W_struct', 'W', 'CovW', 'X', 'CovX', 'Tau_struct');
end

% Load test data
disp('Loading test data..')
data_test = load(['/share/climate/jluttine/testbed/' ...
                  'testbed_vbrfa2011_testdata']);

% Compute predictive log-likelihood
disp('Computing predictive log-likelihood..')
samples = 10;
[M,N] = size(data_test.observations);
Ytest = data_test.observations;
Itest = ~isnan(Ytest);
loglike = nan(samples,1);
for i=1:samples
  
  % Sample a noise level
  tau = gamrnd(Q.W_struct.rho_struct.a_tau, 1/Q.W_struct.rho_struct.b_tau);

  % Sample a reconstruction
  for m=1:M
    W(:,m) = gaussian_rand(Q.W(:,m), 1/tau * Q.CovW(:,:,m));
  end
  for n=1:N
    X(:,n) = gaussian_rand(Q.X(:,n), Q.CovX(:,:,n));
  end
% $$$   W = mvnrnd(Q.W', 1/tau * Q.CovW)';
% $$$   X = mvnrnd(Q.X', Q.CovX)';
  F = W'*X;

  % Log-likelihood
  switch model
   case 1 % GAUSSIAN
    Z = Ytest(Itest) - F(Itest);
    loglike(i) = gaussian_logpdf(tau*(Z'*Z), ...
                                 0, ...
                                 0, ...
                                 -numel(Z)*log(tau), ...
                                 numel(Z));
    
   case 2 % MULTIVARIATE T
    alpha = Q.Tau_struct.a_u;
    beta = Q.Tau_struct.b_u;
    tau = tau * (alpha./beta);
    nu = 2*alpha;
    Z = Ytest - F;
    Z(~Itest) = 0;
    Z2 = tau .* dot(Z,Z,1);
    Nmv = sum(Itest,1)>0;
    loglike(i) = sum(t_logpdf(Z2(Nmv), ...
                              -log(tau(Nmv)), ...
                              nu(Nmv), ...
                              sum(Itest(:,Nmv),1)));
    
   case 3 % INDEPENDENT T
    [I,J] = find(Itest);
    Z2 = tau * (Ytest(Itest) - F(Itest)).^2;
    loglike(i) = sum(t_logpdf(Z2, ...
                              -log(tau), ...
                              Q.Tau_struct.nu(I), ...
                              1));
    
   case 4 % LAPLACE
    loglike(i) = sum(laplace_logpdf(Ytest(Itest), ...
                                    F(Itest), ...
                                    sqrt(tau)));
  end

end
%loglike = loglike / samples;
fprintf('Mean predictive log-density: %.4e\n', mean(loglike));
bias = max(loglike);
fprintf('Log predictive density: %.4e\n', log(mean(exp(loglike-bias))) + bias);

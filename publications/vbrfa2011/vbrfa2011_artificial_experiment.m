
function vbrfa2011_artificial_experiment(seed)

if nargin < 1
  seed = 10;
end
rand('state', seed);
randn('state', seed);

% Observations:
% 
% - When increasing M, Laplace -> independent t
%
% - Variance of corruption makes no big difference between Laplace and
% independent t

% From Bishop, standard deviations 5,4,3,2,1,1,1,1,... :

M = 10;
N = 100;
D = 4;

% Covariance matrix singular values: 

R = orth(randn(M));
eig = ([D:-1:1,zeros(1,M-D)]).^2;
Cov = R*diag(eig)*R';
% $$$ eig = (1+[D:-1:1,zeros(1,M-D)]).^2;
% $$$ Cov = R*diag(eig-1)*R';

mu = randn(M,1);

% Noiseless data
Y = mvnrnd(mu, Cov, N)';

% Noisy data
Yn = Y + 1*randn(M,N);

% Corrupted data
I = rand(M,N) < 0.02;
Yno = Yn + I.*unifrnd(-30,30,M,N);

%
% Construct the model
%

Dh = M-1;

% $$$ tsplot(Yno, 'k')
% $$$ addtsplot(Y, 'b');
% $$$ return

for model=1:4
  
  % It seems that (at least in some cases) it would be better to combine the
  % isotropic noise with X instead W for better and more stable results..?

  noise_iso = noise_module_isotropic(M, 1, ...
                                     'init', struct('a_tau',1, ...
                                                    'b_tau',1));

  % ARD for W
  W_module = factor_module_ard(Dh+1, M, ...
                               'update_alpha', 1, ...
                               'noise_module', noise_iso);

  % IID for X
  mu = [1; zeros(Dh,1)];
  Cov = diag([1e-6; ones(Dh,1)]);
  X_module = factor_module_iid(Dh+1, N, ...
                               'prior', struct('mu', mu, ...
                                               'CovX', Cov), ...
                               'noise_module', []);
                               

  switch model

   case 1 % GAUSSIAN
    noise = noise_module_fixed(M, N, 1);
    id = 'gaussian';

   case 2 % MULTIVARIATE T
    noise = noise_module_multivariate_t(M, N, ...
                                        'update_nu', [1], ...
                                        'init', struct('nu', 0.1));
    id = 'multi-t';

   case 3 % INDEPENDENT T
    noise = noise_module_independent_t(M, N, ...
                                       'update_nu', [1], ...
                                       'nu', 'pooled', ...
                                       'init', struct('nu', 0.1));
    id = 'ind-t';

   case 4 % LAPLACE
    noise = noise_module_laplace(M, N);
    id = 'laplace';
    
   otherwise
    error('Unknown model requested')
    
  end

  %
  % Run the method
  %

  % General function call for running the method
  Q = vbfa(Dh+1, ...
           Yno, ...
           W_module, ...
           X_module, ...
           noise, ...
           'maxiter', 50, ...
           'rotate', [1:10, 10:5:100], ...
           'rotation_maxiter', 20, ...
           'update_x', 1, ...
           'update_w', 1, ...
           'update_noise', 1);

  Yh = Q.W'*Q.X;

  rmse1_train(model) = rmse(Yh, Yno)
  rmse2_test(model) = rmse(Yh(~I), Y(~I))
  rmse3_outliers(model) = rmse(Yh(I), Y(I))

% $$$   tsplot(Yno(1:10,:), 'k')
% $$$   addtsplot(Y(1:10,:), 'b');
% $$$   addtsplot(Yh(1:10,:), 'r');
end


rmse1_train
rmse2_test
rmse3_outliers

%Q.W*diag(Q.rhow)*Q.W' + sum(Q.CovW,3)

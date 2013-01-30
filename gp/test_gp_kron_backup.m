
function test_gp_kron

%randn('state', 100);
%rand('state', 100);

%
% Generate data
%

logtheta1 = log(30);
N1 = 400;
x1 = 1:N1;
K1 = gp_cov_pp(logtheta1, x1, x1);

logtheta2 = log(30);
N2 = 500;
x2 = 1:N2;
K2 = gp_cov_pp(logtheta2, x2, x2);

s = 1e-0;
[X2,X1] = meshgrid(x2,x1);
freq1 = 0.03*2*pi;
freq2 = 0.07*2*pi;
Y_noiseless = lchol(K1) * randn(N1,N2) * lchol(K2)';
Y = Y_noiseless + s*randn(N1,N2);


%
% Get posterior
%

% $$$ % Pseudo code for ideal covariance function construction: (???)
% $$$ covfunc = gp_cov_pp();
% $$$ noisefunc = gp_cov_noise_iso();
% $$$ cov11 = covfunc.set_inputs(x);
% $$$ cov12 = cov11;
% $$$ noise = noisefunc.set_inputs(x);
% $$$ cov22 = @(logtheta) (cov11(logtheta) + noise(logtheta));
% $$$ % or a bit simpler:
% $$$ cov11 = gp_cov_pp(x);
% $$$ cov12 = cov11;
% $$$ noise = gp_cov_noise_iso(x);
% $$$ cov22 = @(logtheta) (cov11(logtheta) + noise(logtheta));
% $$$ cov22 = gp_cov_sum(cov11, noise); % could be a struct?
% $$$ % then by applying some wrapper, you could actually map hyperparameters
% $$$ % arbitrarily to the covariance functions (e.g., use common parameters)
% $$$ %
% $$$ % or if predicting to different locations
% $$$ cov11 = gp_cov_pp(xh);
% $$$ cov12 = gp_cov_pp(xh,x)
% $$$ cov22 = gp_cov_sum( gp_cov_pp(x), gp_cov_noise_iso(x) ); % could be a struct?
% $$$ 
% $$$ % gp_cov_kron(...)
% $$$ % gp_cov_sparse(...)
% $$$ % gp_cov_prod(...)
% $$$ % gp_cov_scale(...)
% $$$ 
% $$$ % maxiter = 100;
% $$$ % [F, theta] = gp_kron_em(Y, covfunc1, covfunc2, maxiter);
% $$$ nsamples = 100;
% $$$ results = gp_kron_mc(Y, covfunc1, covfunc2, 'nsamples', nsamples);
% $$$ 
% $$$ results = gp_mc(y, covfunc_f, covfunc_fy, covfunc_y, logtheta_init, ...
% $$$                 'nsamples', nsamples);


% Missing values
Imv = (rand(N1,N2) < 0.4);
Ymv = Y;
Ymv(Imv) = nan;

covfunc1 = gp_cov_pp_proto(x1);
covfunc2 = gp_cov_pp_proto(x2);
precond1 = [];
precond2 = [];
%precond1 = @(K,s) get_linsolver_ldlchol(K+s*speye(size(K)));
%precond2 = @(K,s) get_linsolver_ldlchol(K+s*speye(size(K)));
F_pcg = 0;
%F_pcg = gp_kron_mc(Ymv,covfunc2,covfunc1,precond2,precond1, logtheta2, ...
%                   logtheta1, s, 'nsamples',1);

theta1 = exp(logtheta1);
theta2 = exp(logtheta2);
F_gibbs = gp_kron_gibbs(Ymv,covfunc2,theta2,1*sqrt(s),covfunc1,theta1, ...
                        1*sqrt(s), 'nsamples',1);

cl = max(abs(Y(:))) + 0.1;
clim = [-cl cl];


Ymv(Imv) = cl;

figure(1)
clf

subplot(3,1,1)
imagesc(Ymv,clim)

subplot(3,1,2)
imagesc(F_pcg,clim)

subplot(3,1,3)
imagesc(F_gibbs,clim)

map_colormap();
cm = colormap();
cm(end,:) = [0.5 0.5 0.5];
colormap(cm);

rmse_Y = rmse(Y_noiseless,Y)
rmse_F = rmse(Y_noiseless,F_pcg)
rmse_Fkron = rmse(Y_noiseless,F_gibbs)



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x1 = gaussian_rand_conditional(x2, K, K1, z1, z2)

if isnumeric(K)
  invK = @(x) K\x;
else
  invK = K;
end
if isnumeric(K1)
  K1 = @(x) K1*x;
end

mu1 = 0;
mu2 = 0;

%x1 = z1;
%z = x1 + z2 - x2 + mu2;
%z = invK(z);                % solves K\z
x1 = mu1 + z1 - K1( invK( z1 + z2 - x2 + mu2 ) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = get_linsolver_ldlchol(K)

LD = ldlchol(K);
M = @(x) linsolve_ldlchol(LD,x);
end

% $$$ function x = mcmc_metropolis_step(x, proposal_rand, target_density)
% $$$ xn = proposal_rand(x);
% $$$ r = target_density(xn) / target_density(x);
% $$$ if rand < r
% $$$   x = xn;
% $$$ end
% $$$ 
% $$$ end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = gp_kron_gibbs(Y, ...
                           covfunc1, theta1, s1, ...
                           covfunc2, theta2, s2, ...
                           varargin)

options = struct( ...
    'nsamples',    1);

Imv = isnan(Y);
[N2,N1] = size(Y);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Initialize missing values
sz = size(Y)
Y(Imv) = 0;

ngibbs = 100

disp('Using joint kronecker cov (gibbs)..')

sample_y = get_sampler_y(covfunc1,theta1,s1,covfunc2,theta2,s2)

thetas = zeros(numel([theta1(:);s1(:);theta2(:);s2(:)]), ngibbs);
tic
for n=1:ngibbs
  
  fprintf('Iteration step %d\n', n);
  
  %
  % Sample global stuff
  %
  
  % Sample A (and its hyperparameters)
  
  % Sample S (and its hyperparameters)
  
  %
  % Sample short-scale stuff
  %
  
  % Sample F (and its hyperparameters)
  
  [Y,F,thetas(:,n)] = sample_y(Y,Imv);
  
  %
  % Sample missing data
  %
  
  % Sample Y
  %Y(~Iobs) = F_noisy(~Iobs) + s*randn(sum(~Iobs(:)),1);
  
  
end
toc

figure(2)
semilogy(thetas');

% F = Y;

% If you want the noisy reconstructions:
% F = F_noisy;

end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function y_sampler_func = get_sampler_y(covfunc1, theta1, s1, covfunc2, ...
                                                    theta2, s2)
  
  % gp_
  
  % Takes:
  % - a function handle covfunc which returns the covmat-struct
  % - a function handle get_sampler which returns sampler function handle
  % - a function handle sample_hyperparameters?
  % - a function handle likelihood?
  
  % this system tries to fulfill the requirements of the Gibbs sampler
  % for GP + FA model - don't make it more general.. the main purpose is
  % to model the residuals and sample missing observations
  
  % Options:
  %
  % - proposal distribution for theta
  
  n_theta1 = numel(theta1);
  n_theta2 = numel(theta2);
  theta = [theta1(:); s1; theta2(:); s2];
  covmatrix = struct();
  sample_variables = [];

  n = 0;
  Ymean = 0;
  
  % Initialize stuff with the given logtheta initialization
  covmatrix = covfunc(theta);
  sample_variables = get_sampler();
  theta
  
  y_sampler_func = @y_sampler;
  
  % returns y-sampler and some solver for FA model?
  covmatrix_linsolver = [];
  
    function [Y, F_noiseless, thetas] = y_sampler(Y,I)
    sample_hyperparameters(Y);
    [F, F_noiseless] = sample_variables(Y);
    % Reconstruct missing values
    if nargin < 2
      Y = F + covmatrix.s1*covmatrix.s2*randn(size(F));
    else
      Y(I) = F(I) + covmatrix.s1*covmatrix.s2*randn(size(F(I)));
    end
    thetas = theta;
    
    Ymean = (n*Ymean + F) / (n+1);
    Ymean = (n*Ymean + F_noiseless) / (n+1);
    n = n+1;
    F_noiseless = Ymean;
    
    %F_noiseless = F;
    end

    function sample_hyperparameters(Y)
    % Sample new hyperparameters using Metropolis
    theta_prop = proprand(theta);
    covmatrix_prop = covfunc(theta_prop);
    if log(rand) < (targetdist(covmatrix_prop,Y) - targetdist(covmatrix,Y))
      covmatrix = covmatrix_prop;
      theta = theta_prop;
      sample_variables = get_sampler();
      disp('accepted')
    else
      disp('rejected')
    end
    end

    function theta = proprand(theta)
    theta = exp(log(theta) + 0.001*randn(size(theta)));
    end

    function covmat = covfunc(theta)
    theta1 = theta(1:n_theta1);
    covmat.K1 = covfunc1(log(theta1));
    covmat.s1 = theta((n_theta1+1));
    covmat.LD_Cov1 = ldlchol(covmat.K1 + covmat.s1^2*speye(size(covmat.K1)));
    
    theta2 = theta((n_theta1+2):(n_theta1+1+n_theta2));
    covmat.K2 = covfunc2(log(theta2));
    covmat.s2 = theta(end);
    covmat.LD_Cov2 = ldlchol(covmat.K2 + covmat.s2^2*speye(size(covmat.K2)));
    end
  
    function logp = targetdist(covmat, Y)
    % gaussian_logpdf_kron_ldlchol
    logp = gaussian_logpdf(...
        Y(:)'*linsolve_kron_ldlchol(covmat.LD_Cov1,covmat.LD_Cov2,Y,'vector'), ...
        0, ...
        0, ...
        logdet_kron_ldlchol(covmat.LD_Cov1,covmat.LD_Cov2), ...
        numel(Y));
    end

    function sampler_func = get_sampler()
    L1 = lchol(covmatrix.K1);
    L2 = lchol(covmatrix.K2);
    N1 = size(covmatrix.K1,1);
    N2 = size(covmatrix.K2,1);
    s1 = covmatrix.s1;
    s2 = covmatrix.s2;
    LD_Cov1 = covmatrix.LD_Cov1;
    LD_Cov2 = covmatrix.LD_Cov2;
    K1 = covmatrix.K1;
    K2 = covmatrix.K2;
    Covfun = @(Y) linsolve_kron_ldlchol(LD_Cov1,LD_Cov2,Y);
    Kfun_noisy = @(X) (kronprod(K1,K2,X) + ...
                       kronprod(K1,s2^2,X) + ...
                       kronprod(s1^2,K2,X));
    Kfun = @(X) kronprod(K1,K2,X);
    sampler_func = @sampler;
           
      function [F_noisy, F_noiseless] = sampler(Y)
      % function [F_noisy, F] = sampler(Y)
      z_k1k2 = kronprod(L1,L2,randn(N2,N1));
      z_k1s_sk2 = kronprod(L1,s2,randn(N2,N1)) + ...
          kronprod(s1,L2,randn(N2,N1));
      z_s = s1*s2*randn(N2,N1);
      F_noiseless = gaussian_rand_conditional(Y, Covfun, Kfun, z_k1k2, ...
                                              z_k1s_sk2+z_s);
      F_noisy = gaussian_rand_conditional(Y, Covfun, Kfun_noisy, z_k1k2+z_k1s_sk2, ...
                                          z_s);
      end
    end

  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = gp_kron_mc(Y, covfunc1, covfunc2, precond1, precond2, ...
                        hyperparams1, hyperparams2, s, varargin)

options = struct( ...
    'nsamples',    100);


Iobs = ~isnan(Y);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% hyperparameters
[K1,L1] = covfunc1(hyperparams1);
[K2,L2] = covfunc2(hyperparams2);
s;

if false
  tic
    % This gets out of memory even for > 1000x1000 K1/K2..
    disp('Cholesky approach..')
    nnz1 = nnz(K1)
    nnz2 = nnz(K2)
    K = kron(K1,K2);
    nnz_joint = nnz(K)
    LD = ldlchol(K);
  toc
end

for n=1:options.nsamples

  %
  % Sample F
  %
  
  % Get preconditioner
  if ~isempty(precond1)
    M1 = precond1(K1,s);
  else
    M1 = [];
  end
  if ~isempty(precond2)
    M2 = precond2(K2,s);
  else
    M2 = [];
  end
  
  % Draw the sample
  tic
    disp('PCG sampling')
    F = gaussian_kron_rand(L1, L2, K1, K2, Y, s, s^2, M1, M2, Iobs);
  toc
  
  %
  % Sample hyperparameters
  %
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COV(F) = KRON(K1,K2)
% COV(Y|F) = K_NOISE
% K1 = L1*L1'
% K2 = L2*L2'
% K_NOISE = L_NOISE*L_NOISE'
% INV(KRON(M1,M2)) approximates INV(KRON(K1,K2)+K_NOISE)
% K_NOISE, L_NOISE, M1 and M2 can be function handles.
function F = gaussian_kron_rand(L1, L2, K1, K2, Y, L_noise, K_noise, M1, ...
                                M2, mask)

% Rename to, e.g., gaussian_kron_rand_conditional

% TODO:
% - missing values

% In order to make this more general, I might have to remove L1, L2 and
% L_noise and use kron(L1,L2)*randn(...) and L_noise*randn(...) instead,
% so one can utilise the structure of the matrices better.

N1 = size(L1,2);
N2 = size(L2,2);

% $$$ U1 = L1';
% $$$ U2 = L2';

if nargin == 2

  % Sample from the prior
  F = kronprod(L1,L2,randn(N2,N1));

else

  % Sample from the posterior
  
  I = speye(N1*N2);
  P = I(mask(:),:);
  
  if isnumeric(L_noise)
    L_noise = @(x) L_noise*x;
  end
  z_noise = L_noise(randn(N1*N2,1));
  
  if isnumeric(K_noise)
    K_noise = @(x) K_noise*x;
  end

  z_f = kronprod(L1,L2,randn(N2,N1),'vector');
% $$$   L_f = @(x) kronprod(L1,L2,x,'vector');
  K_f = @(x) kronprod(K1,K2,reshape(x,N2,N1),'vector');
  K_y = @(x) (K_f(x) + K_noise(x));
  if ~isempty(M1) && ~isempty(M2)
    M_y = @(x) linsolve_kron(M1, M2, reshape(x,N2,N1), 'vector');
  elseif ~isempty(M2)
    error('not yet implemented')
  elseif ~isempty(M1)
    error('not yet implemented')
  else
    M_y = [];
  end
  
  F = gaussian_rand_pcg(Y(:), K_y, K_f, z_f, z_noise, M_y, P, 'maxiter', ...
                        100, 'tol', 1e-3);
% $$$   F = gaussian_rand_pcg(Y(:), K_y, K_f, L_f, L_noise, M_y);
  F = reshape(F,N2,N1);
  
end

end



















%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ function results = gp_mc(y, covfunc, logtheta0, varargin)
% $$$ 
% $$$ options = struct( ...
% $$$     'nsamples',    100);
% $$$ 
% $$$ % Parse arguments
% $$$ [options, errmsg] = argparse( options, varargin{:} );
% $$$ error(errmsg);
% $$$ 
% $$$ covmatrix = covfunc.set_hyperparameters(logtheta); % a struct
% $$$ 
% $$$ for n=1:options.nsamples
% $$$   
% $$$   % Sample function values of interest
% $$$   % [Cov_f, Cov_fy, Cov_y] = covfunc.get_covariance(logtheta);
% $$$   % covfunc.
% $$$   %f = gp_draw_sample(Cov_f, y, Cov_y, Cov_fy);
% $$$   f = options.sample_f(covmatrix, y);
% $$$   
% $$$   % Sample hyperparameters
% $$$   hyperparams = options.sample_hyperparameters(covfunc, f, y, hyperparams);
% $$$   covmatrix = covfunc.set_hyperparameters(hyperparams);
% $$$   
% $$$ end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ function f = gp_draw_sample(Cov_f, y, Cov_y, Cov_fy)
% $$$ % $$$ C_f = Cov_f.get_covariance_matrix();
% $$$ % $$$ if nargin >= 2
% $$$ % $$$   L_y = Cov_y.get_covariance_matrix();
% $$$ % $$$   C_fy = Cov_fy.get_covariance_matrix();
% $$$ % $$$   L = chol(C_f - C_fy
% $$$ 
% $$$ L_y = chol(Cov_y,'lower');
% $$$ K = linsolve_tril(L_y,Cov_fy');
% $$$ L = chol(Cov_f - K'*K);
% $$$ mu = K'*linsolve_tril(L_y,y);
% $$$ f = gaussian_rand(mu,L);
% $$$   
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ % E(X) = MU
% $$$ % COV(X) = L*L' (L does not need to be square)
% $$$ function x = gaussian_rand(mu, L)
% $$$ x = mu + L*randn([size(L,2),1]);
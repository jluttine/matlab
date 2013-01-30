
function gp_init_kron(covfun1, covfun2)

% low-level approach..?

cov1 = gp_cov_pp();
cov2 = gp_cov_noise_iso();
cov = gp_cov_sum(cov1,cov2);

gp_moments()

function [mu, Cov] = gp_moments(y,Ky,Kxy,Kx)

% mu = K_xy * (K_y \ y);
% Cov = K_x - K_xy * inv(K_y) * K_yx

if isnumeric(Kxy)
  Kxy = @(x) Kxy*x;
end
if isnumeric(Ky)
  Ky = @(x) linsolve_cov(Ky,x);
end

mu = Kxy( Ky( y ) );
if nargout >= 2
  if isvector(Kx)
    % Evaluate variances
    Cov = Kx - 
  else
    % Evaluate covariance
    Kxy = Kxy(speye(size(Ky1)));
    Cov = Kx - Kxy
  end
end

end

% regular gp:

% set up covariance function
%
% hyperparameter stuff should be removed.. GP is always conditional on
% the given hyperparams?
cov1 = gp_cov_pp(pdf_gamma(2,3), pdf_gamma(0.5,0.5), 'jitter', 1e-6); % struct
gp1 = gp_init(cov1);
cov2 = gp_cov_noise_iso(pdf_gamma(1,1), 'jitter', 0); % struct
gp2 = gp_init(cov2);
% cov = gp_cov_sum_sparse(cov1, cov2);

gp_y = gp_sum(gp1,gp2); % ?

gp_y.set_data_inputs(x);


gp = gp_init(cov)

gp.set_data_inputs(x)

gp.set_hyperparameters(logtheta, ..)
sampler = gp.get_hyperparameter_sampler(..)

gp.get_posterior_moments(logtheta,..)
sampler = gp.get_posterior_sampler(logtheta, ..)
f = sampler(y, ..)
gp.get_hyperparameter_optimizer(..)




% or:
%gp = gp_init(prior_covariance, likelihood_covariance);




%%%%

gp.set_data_inputs(..)

gp.optimize_hyperparameters(..) % ?
gp.draw_hyperparameters(..) % 

gp.evaluate_prior_covariance(hyperparameters)
gp.evaluate_likelihood_covariance(noise_covariance)
gp.draw_sample(y); % without data draws from prior

%
% CONSTRUCT THE GP
%

% Create covariance structs
cov1 = gp_init_cov_pp(); % struct
cov1.params.prior = pdf_product(pdf_gamma(2,3), pdf_gamma(0.5,0.5)); % struct
cov1.jitter = 1e-6;
% functions:
% cov1.data_inputs(x) - evaluates, e.g., distance matrix
% cov1.evaluate
% cov1.cholesky()
% cov1.linsolve()
% cov1.multiply()
cov2 = gp_init_cov_pp(); % struct
cov2.params.prior = pdf_product(pdf_gamma(2,3), pdf_gamma(0.5,0.5)); % struct
cov2.jitter = 1e-6;

% Create Kronecker product covariance struct
cov = gp_cov_kron(cov1, cov2);

% Also, there could be, e.g., mean function struct?

gp = gp_init(cov); % ????

%
% ASSIGN DATA
%

gp = gp_set_data_inputs(gp, {x1,x2}); % makes it possible to pre-evaluate,
                                      % e.g., distance matrix

gp_set_hyperparameters() % ??

gp_evaluate_prior_covariance(); % evaluate K and possibly L_K ?

gp_observations(Y,s^2);

% when to evaluate preconditioner?

%
% DO INFERENCE
%

gp_sample_cg(); %% ??

f = gp_sample();

% you could also have functions such as gp_mean, gp_variance for VB
% inference..?i






N1 = [];
N2 = [];

K1 = [];
K2 = [];

LD1 = [];
LD2 = [];

L1 = [];
L2 = [];

N1 = 500;
x1 = 1:N1;
K1 = gp_cov_pp(log(10), x1, x1);

N2 = 500;
x2 = 1:N2;
K2 = gp_cov_pp(log(10), x2, x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function set_parameters(logt1,logt2,x1,x2)
  logtheta1 = logt1;
  X1 = x1;
  N1 = size(x1,2);
  
  logtheta2 = logt2;
  X2 = x2;
  N2 = size(x2,2);
  
  evaluate_covariance();
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function rndfunc = get_posterior_sampling_function(Y,noise_std)
  
  s = noise_std;

  Afun = @(x) (kronprod(K1,K2,x,'vector')+s^2*x);
  Mfun = get_preconditioner();
  
  rndfunc = @() draw_sample(Y,Afun,Mfun);
  
  end

  function F = draw_sample()
  
  Z0 = randn(N1,N2);
  Z = randn(N1,N2);
  F0 = kronprod(L2,L1,Z0);
  G = F0 + s*Z - Y;
  X = reshape( pcg(Afun, G(:), 1e-4, 20, Mfun), N1, N2 );
  F = F0 - kronprod(K1,K2,X);
  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function precondfun = get_preconditioner()
  LD_C1 = ldlchol(K1+s*speye(N1));
  LD_C2 = ldlchol(K1+s*speye(N1));
  precondfun = @(x) kronsolve_for_cg2(LD_C1, LD_C2, x);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function evaluate_covariance()
  K1 = covfun1(logtheta1, X1, X1);
  LD1 = ldlchol(K1);
  [L,D] = ldlsplit(LD1);
  L1 = L*sqrt(D);
  
  K2 = covfun2(logtheta2, X2, X2);
  LD2 = ldlchol(K2);
  [L,D] = ldlsplit(LD2);
  L2 = L*sqrt(D);
  end
  
end
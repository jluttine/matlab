
% [covmatrix] = covfunc(theta)
% [covmatrix, theta] = rand_theta(theta, covmatrix, Y, covfunc,
%                                 logposterior)
% [rand_y] = get_randfunc_y(covmatrix)
% Y = rand_y(covmatrix, Y, I)
% [y, dy] = logprior_theta(theta)
% [y, dy] = loglikelihood_theta(Y, covmatrix)

% COVMATRIX:
% L1/LD1
% L2/LD2
% linsolve1
% linsolve2
% K1
% K2
% logdet1
% logdet2

function res = test_gp_get_randfunc()

%
% DATA
%

%randn('state', 10)
%rand('state', 10)

% Inputs
N1 = 14;
x1 = 1:N1;
N2 = 20;
x2 = 1:N2;

% Covariance models
d1 = sqrt(sq_dist(x1(1), x1));
covfunc1 = gp_cov_pp(d1,1);
covfunc1 = gp_cov_scale(covfunc1);
covfunc1 = gp_cov_toeplitz(covfunc1);
theta1 = [1; 20];
K1 = covfunc1(theta1);
d2 = sqrt(sq_dist(x2(1), x2));
covfunc2 = gp_cov_pp(d2,1);
covfunc2 = gp_cov_scale(covfunc2);
covfunc2 = gp_cov_toeplitz(covfunc2);
theta2 = [1; 20];
K2 = covfunc2(theta2);

% $$$ % Plot gradients
% $$$ [K1, dK1] = covfunc1(theta1);
% $$$ figure
% $$$ subplot(3,1,1)
% $$$ imagesc(K1)
% $$$ subplot(3,1,2)
% $$$ imagesc(dK1{1})
% $$$ subplot(3,1,3)
% $$$ imagesc(dK1{2})
% $$$ return


% Generate noisy data
s = 10e-1;
Y = lchol(K2) * randn(N2,N1) * lchol(K1)';
save('Y_noiseless', 'Y');
Y = Y + s*randn(N2,N1);
save('Y_noisy', 'Y');

% Missing values
pmv = 0.2;
Imv = randperm(N2*N1);
Imv = Imv(1:floor(pmv*(N2*N1)));
%Imv = (rand(N2,N1) < 0.2);
Y(Imv) = nan;

%
% INFERENCE
%

% ??
% covsum = gp_cov_sum(@gp_cov_pp, @gp_cov_pp);
% covfunc = covsum({D1,d}, {D2,d})
% K = covfunc({30}, {40})
% covfunc = gp_cov_mapper_simple(covfunc);
% K = covfunc([30 40]);
% ?? (K might not always be the covariance matrix but something else?)
%
% covscale = @(D) gp_cov_scale(covfunc(D));
% covsum = @(D) gp_cov_sum( covfunc_a(D), gp_cov_scale(covfunc_b(D)) );


% covfunc1 = @(theta) gp_cov_pp(theta, D1, d);
% covfunc2 = @(theta) gp_cov_pp(theta, D2, d);
% covfuncsum = gp_cov_sum(covfunc1, covfunc2);
% covfuncprod = gp_cov_prod(covfunc1, covfunc2);
% covfuncscale = gp_cov_scale(covfunc); ??
% covfunc = gp_cov_fix(covfunc, [2], [30]); % fixes theta(2)=30
% covfunckron = ??;

% Gibbs sampling for Y(missing) and covariance parameters
N_samples = 100;

switch 'kron2'
 case 'kron1'

  %logprior_theta = @(theta) sum(lognormal_logpdf(theta, 0, 1));
  logprior_theta = @(theta) sum(gamma_logpdf(theta, 1e-3, 1e-3));
  

  % Initial guess for covariance parameters
  theta_init = [1         ... % signal magnitude
                theta1(2) ... % length scale
                sqrt(s)   ... % noise magnitude
                1         ... % signal magnitude
                theta2(2) ... % length scale
                sqrt(s)]';    % noise magnitude
  
  samplefunc = get_sample_function2(numel(theta_init), N_samples);

  func1 = @get_logposterior;
  func2 = @get_randfunc;
  func3 = @covstructfunc;
  [get_logpdf, get_rand_model, func_theta] = gp_init_kron1(covfunc1, ...
                                                    @get_covsolver_ldlchol,    ...
                                                    @(K) chol(K,'lower'),      ...
                                                    covfunc2,                  ...
                                                    @get_covsolver_ldlchol,    ...
                                                    @(K) chol(K,'lower'),      ...
                                                    logprior_theta,            ...
                                                    samplefunc);

 case 'kron2'
   
  %logprior_theta = @(theta) sum(gamma_logpdf(theta, 1e-3, 1e-3));
  logprior_theta = @(theta) sumout(@() gamma_logpdf(theta, 1e-3, 1e-3));
% $$$   logprior_theta = @(theta) sumout(@() gamma_logpdf(theta(1), 1e-3, 1e-3), ...
% $$$                                    @() gamma_logpdf(theta(2), 1e-3, 1e-3), ...
% $$$                                    @() gamma_logpdf(theta(3), 1e-3, 1e-3), ...
% $$$                                    @() gamma_logpdf(theta(4), 1e-3, 1e-3), ...
% $$$                                    @() gamma_logpdf(theta(5), 1e-3, 1e-3));
  
% $$$   logprior_theta = joint_logpdf({@gamma_logpdf, 1e-3, 1e-3}, ...
% $$$                                 {@gamma_logpdf, 1e-3, 1e-3}, ...
% $$$                                 {@gamma_logpdf, 1e-3, 1e-3}, ...
% $$$                                 {@gamma_logpdf, 1e-3, 1e-3}, ...
% $$$                                 {@gamma_logpdf, 1e-3, 1e-3});
  

  % Initial guess for covariance parameters
  theta_init = [1         ... % signal magnitude
                theta1(2) ... % length scale
                1         ... % signal magnitude
                theta2(2) ... % length scale
                s]';          % noise magnitude
  
  samplefunc = get_sample_function2(numel(theta_init), N_samples);

  [get_logpdf, get_rand_model, func_theta] = gp_init_kron2(covfunc1, ...
                                                    @get_covsolver_ldlchol, ...
                                                    @(K) chol(K,'lower'), ...
                                                    covfunc2, ...
                                                    @get_covsolver_ldlchol, ...
                                                    @(K) chol(K,'lower'), ...
                                                    logprior_theta, ...
                                                    samplefunc);

end


% Posterior sampling for covariance parameters
switch 'slice'
  
 case 'mh'
  % Metropolis-Hastings
  stepsize = 1e-2;
  rand_proposal = @(theta) lognormal_rand(log(theta), stepsize);
  logpdf_proposal = @(theta1, theta2) ...
      sum(lognormal_logpdf(theta1,      ...
                           log(theta2), ...
                           stepsize));
  rand_theta = mcmc_init_metropolishastings(theta_init,      ...
                                            get_logpdf,    ...
                                            rand_proposal,   ...
                                            logpdf_proposal, ...
                                            func_theta);
  
 case 'hmc'
  % Hamiltonian/Hybrid Monte Carlo
  
 case 'slice'
  % Slice sampling (do sampling in log-scale)
  %
  % NOTE: Slice sampling can be very bad if initialized at very low
  % density values. Sometimes it might be surprising that a "good"
  % initial guess happens to have very low density so the algorithm
  % won't work in practice.
  rand_theta = mcmc_init_slicesampling(log(theta_init), ...
                                       get_logpdf, ...
                                       'f', func_theta, ...
                                       'y', @(logtheta) exp(logtheta), ...
                                       'logdy', @(logtheta) sum(logtheta), ...
                                       'dy', @(logtheta) exp(logtheta), ...
                                       'ddy', @(logtheta) exp(logtheta), ...
                                       'type', 'reflective-outside');
  
  function [get_logpdf, func_logtheta] = theta2logtheta(get_logpdf, func_theta)
  y = @(logtheta) exp(logtheta);
  logdy = @(logtheta) sum(logtheta);
    function [f_logtheta, logpdf_logtheta] = get_logpdf(logtheta)
    theta = y(logtheta);
    f_theta = func_theta(theta);
    logpdf_theta = get_logpdf(theta, f_theta);
    logdy_logtheta = logdy(logtheta);
    logpdf_logtheta = @(varargin) (logpdf_theta(varargin{:}) + logdy_logtheta);
    end
  end
  
 case 'ml'
  % Maximum likelihood
  rand_theta = init_maximumlikelihood();
  
 case default
  error('Unknown sampling scheme for theta');

end


rand_model = get_rand_model(rand_theta, theta_init);

% Gibbs sampling
tic
gibbs(Y,Imv, N_samples, rand_model);
toc

%
% ANALYSIS
%

% Check the results
res = samplefunc()

cl = max(abs(Y(:))) + 0.1;
clim = [-cl cl];

Y(Imv) = cl;

figure(1)
clf

subplot(2,2,1)
imagesc(Y,clim)

load('Y_noiseless')

subplot(2,2,3)
imagesc(Y,clim)

subplot(2,2,4)
imagesc(res.F,clim)

map_colormap();
cm = colormap();
cm(end,:) = [0.5 0.5 0.5];
colormap(cm);

rmse_F = rmse(Y,res.F)

subplot(2,2,2)
semilogy(res.theta');

end

% $$$   function [f, df] = sumout(varargin)
% $$$   f = 0;
% $$$   for n=1:length(varargin)
% $$$     if nargout <= 1
% $$$       f = f + varargin{n}();
% $$$     else
% $$$       [fn, dfn] = varargin{n}();
% $$$       f = f + fn;
% $$$       df = [df; dfn];
% $$$     end
% $$$   end
% $$$   end

function [f, df] = sumout(g)
if nargout <= 1
  f = sum(g());
else
  [f, df] = g();
  f = sum(f);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function covfunc_toeplitz = gp_cov_toeplitz(covfunc)
covfunc_toeplitz = @cov_toeplitz;
  function varargout = cov_toeplitz(theta)
  varargout = cell(nargout,1);
  if nargin == 0
    if nargout == 1
      varargout{1} = covfunc();
    else
      out = cell(max(nargout,3),1);
      [out{:}] = covfunc();
      out{2} = max(out{2}, out{3});
      out{3} = out{2};
      varargout(:) = out(1:nargout);
    end
    return
  end
  if nargout == 1
    k = covfunc(theta);
    varargout{1} = toeplitz(k);
  elseif nargout == 2
    [k, dk] = covfunc(theta);
    varargout{1} = toeplitz(k);
    varargout{2} = cell(numel(dk),1);
    for i=1:numel(dk)
      varargout{2}{i} = toeplitz(dk{i});
    end
  end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func = gp_cov_scale(covfunc, dcovfunc)
func = @cov;
  function varargout = cov(theta)
  out = cell(nargout,1);
  varargout = cell(nargout,1);
  if nargin == 0
    [out{:}] = covfunc();
    varargout{1} = 1 + out{1};
    varargout(2:end) = out(2:end);
    return
  end
  [out{:}] = covfunc(theta(2:end));
  varargout{1} = theta(1)^2 * out{1};
  if nargout >= 2
    varargout{2} = cell(numel(theta),1);
    varargout{2}{1} = 2*theta(1)*out{1};
    for n=2:numel(theta)
      % Scale the gradients
      varargout{2}{n} = theta(1)^2 * out{2}{n-1};
    end
  end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samplefunc = get_sample_function2(D_theta, N)
results.Y = [];
results.F = 0;
results.theta = zeros(D_theta, N);
samplefunc = @process_sample;
n = 1;
  function res = process_sample(Y, F, theta)
  if nargin >= 1
    results.Y = Y;
    switch 'sample'

     case 'mean'
      results.F = (F + (n-1)*results.F) / n;
      
     case 'sample'
      results.F = F;
      
     case default
      error('Unknown storing');
      
    end
    results.theta(:,n) = theta(:);
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
  
  Y = rand_model(Y,Imv);
  fprintf('Iteration step %d done.\n', n)
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [func1, func2, func3] = gp_init_kron1(covfunc1, ...
                                               get_solver1,     ...
                                               get_squareroot1, ...
                                               covfunc2,        ...
                                               get_solver2,     ...
                                               get_squareroot2, ...
                                               logprior,        ...
                                               samplefunc)
% additive noise to both kron-components:
% kron(K1+s1^2, K2+s2^2)

[n_theta1, N1] = covfunc1();
[n_theta2, N2] = covfunc2();

% Sampler for theta
get_loglikelihood = @get_loglikelihood_marginalized;

func1 = @get_logposterior;
func2 = @get_randfunc;
func3 = @covstructfunc;

  function func = get_logposterior(theta, covstruct, covstruct_grad)
  %lprior = 0;
  if nargin <= 2
    loglikelihood = get_loglikelihood(covstruct);
    lprior = logprior(theta);
  else nargin >= 3
    loglikelihood = get_loglikelihood(covstruct, covstruct_grad);
    [lprior, dlprior] = logprior(theta);
  end
  func = @logposterior;
    function [lp, dlp] = logposterior(Y,varargin)
    lp = lprior + loglikelihood(Y,varargin{:});
    end
  end
  
  function func = get_randfunc(rand_theta, theta_init)
  % Initialize latent variables
  F = [];%zeros(N2,N1);
  covstruct = covstructfunc(theta_init);
  get_rand_y = @get_rand_y_gibbs;
  %get_rand_y = @get_rand_y_pcg;
  rand_y = get_rand_y_pcg(covstruct);

  func = @randfunc;
    function Y = randfunc(Y,I)
  
    if isempty(F)
      % Sample Y
      [Y, F] = rand_y(Y,I);
    end

    % Sample hyperparameters
    [theta, covstruct_new] = rand_theta(Y);
    if ~isequal(covstruct, covstruct_new)
      covstruct = covstruct_new;
      rand_y = get_rand_y_gibbs(covstruct);
    end
    
    % Sample Y
    [Y, F] = rand_y(Y,I);
    
    % Process samples??
    if ~isempty(samplefunc)
      samplefunc(Y, F, theta);
    end
    end
  end

  function [covstruct, covstruct_grad] = covstructfunc(theta, varargin)
  
  theta1 = theta(1:n_theta1);
  covstruct.K1 = covfunc1(theta1);
  covstruct.s1 = theta(n_theta1+1);
%  [linsolve1, logdet1] = get_solver1(covstruct.K1);
  [linsolve1, logdet1] = get_solver1(covstruct.K1, covstruct.s1^2);
  % multiply1 = get_multiplier1(K1);

  theta2 = theta((n_theta1+2):(end-1));
  covstruct.K2 = covfunc2(theta2);
  covstruct.s2 = theta(end);
  %[linsolve2, logdet2] = get_solver2(covstruct.K2);
  [linsolve2, logdet2] = get_solver2(covstruct.K2, covstruct.s2^2);
  % multiply2 = get_multiplier2(K2);

  

  covstruct.multiply_f = @(X) (kronprod(covstruct.K1,covstruct.K2,X));

  covstruct.multiply_noiseless = @(X) (kronprod(covstruct.K1,covstruct.K2,X) + ...
                                       kronprod(covstruct.K1,covstruct.s2^2,X) + ...
                                       kronprod(covstruct.s1^2,covstruct.K2,X));

  covstruct.multiply = @(X) (covstruct.multiply_noiseless(X) + ...
                             covstruct.s1^2*covstruct.s2^2*X);
  
  covstruct.linsolve = @(X) linsolve_kron(linsolve1, linsolve2, X);
  
  covstruct.logdet = @() N2*logdet1() + N1*logdet2();

  end

  function func = get_rand_y_gibbs(covstruct)

  L1 = get_squareroot1(covstruct.K1);
  L2 = get_squareroot2(covstruct.K2);
  func = @rand_y_gibbs;

    function [Y, F] = rand_y_gibbs(Y, Imv)
    % Sample latent function values
    F0 = gaussian_rand_kron(L1,L2);
    F0_1 = gaussian_rand_kron(L1,covstruct.s2*speye(size(L2)));
    F0_2 = gaussian_rand_kron(covstruct.s1*speye(size(L1)),L2);
    Y0 = F0 + F0_1 + F0_2 + covstruct.s1*covstruct.s2*randn(N2,N1);
    Z = covstruct.linsolve( Y0 - Y );
    F = F0 + F0_1 + F0_2 - covstruct.multiply_noiseless(Z);
    % Sample missing values
    Y(Imv) = F(Imv) + covstruct.s1*covstruct.s2*randn(size(Y(Imv)));
    % Interesting function values
    F = F0 - covstruct.multiply_f(Z);
    % You could evaluate:
    % rand(F), mean(F), mean(Y), 
    end
  
  end

  function func = get_rand_y_pcg(covstruct)
  L1 = get_squareroot1(covstruct.K1);
  L2 = get_squareroot2(covstruct.K2);
  func = @rand_y_pcg;

    function [Y,F] = rand_y_pcg(Y, Imv)
    S = speye(numel(Y));
    S(Imv(:),:) = [];
    multiply = @(x) reshape(covstruct.multiply(reshape(x,size(Y))), ...
                            [numel(Y),1]);
    multiply_noiseless = @(x) reshape(covstruct.multiply_noiseless(reshape(x,size(Y))), ...
                              [numel(Y),1]);
    F0 = gaussian_rand_kron(L1,L2);
    F0_1 = gaussian_rand_kron(L1,covstruct.s2*speye(size(L2)));
    F0_2 = gaussian_rand_kron(covstruct.s1*speye(size(L1)),L2);
    %Y0 = F0 + F0_1 + F0_2 + covstruct.s1*covstruct.s2*randn(N2,N1);
    %Z_f = gaussian_rand_kron(L1,L2);
    Z_noise = covstruct.s1*covstruct.s2 * randn(numel(Y),1);
    F = gaussian_rand_pcg(Y(:), ...
                          multiply, ...
                          multiply_noiseless, ...
                          F0(:) + F0_1(:) + F0_2(:), ...
                          Z_noise(:), ...
                          [], ...
                          S, ...
                          'maxiter', 1e3, ...
                          'tol', 1e-5);
    F = reshape(F, size(Y));
    Y(Imv) = F(Imv) + covstruct.s1*covstruct.s2*randn(size(Y(Imv)));
    multiply_f = @(x) reshape(covstruct.multiply_f(reshape(x,size(Y))), ...
                              [numel(Y),1]);
    F = gaussian_rand_pcg(Y(:), ...
                          multiply, ...
                          multiply_f, ...
                          F0(:) + F0_1(:) + F0_2(:), ...
                          Z_noise(:), ...
                          [], ...
                          S, ...
                          'maxiter', 1e3, ...
                          'tol', 1e-5);
    end
  end

  function func = get_loglikelihood_marginalized(covstruct, covstruct_grad)
  % marginalized
  % other options: whitened, conditioned
  ldet = covstruct.logdet();
  func = @loglikelihood_marginalized;
    function [x, dx] = loglikelihood_marginalized(Y)
    Z = covstruct.linsolve(Y);
    x = gaussian_logpdf(...
        Y(:)'*Z(:), ...
        0, ...
        0, ...
        ldet, ...
        numel(Y));
    if nargout >= 2
      for n=1:N
        V = gradient{n}(Z);
        dx(n) = 0.5 * (V(:)'*Z(:)) - ...
                0.5 * trace_linsolve_gradient{n}();
        
      end
      error('hmm hmm');
    end
    end
  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [func1, func2, func3] = gp_init_kron2(covfunc1, ...
                                               get_solver1, ...
                                               get_squareroot1, ...
                                               covfunc2, ...
                                               get_solver2, ...
                                               get_squareroot2, ...
                                               logprior, ...
                                               samplefunc)

% additive noise to both kron-components:
% kron(K1+s1^2, K2+s2^2)


[n_theta1, N1, M1] = covfunc1();
[n_theta2, N2, M2] = covfunc2();

% Posterior density function
get_loglikelihood = @get_loglikelihood_conditioned;
get_logpdf = @get_logposterior;

func1 = @get_logposterior;
func2 = @get_randfunc;
func3 = @covstructfunc;


  function func = get_logposterior(theta, covstruct, covstruct_grad)
  if nargin <= 2
    loglikelihood = get_loglikelihood(covstruct);
    lprior = logprior(theta);
  else % nargin >= 3
    loglikelihood = get_loglikelihood(covstruct, covstruct_grad);
    [lprior, dlprior] = logprior(theta);
  end
  func = @logposterior;
    function [lp, dlp] = logposterior(Y,varargin)
    if nargout <= 1
      lp = lprior + loglikelihood(Y,varargin{:});
    elseif nargout == 2
      [llike,dllike] = loglikelihood(Y,varargin{:});
      lp = llike + lprior;
      dlp = dllike + dlprior;
    end
    end
  end
  
  function func = get_randfunc(rand_theta, theta_init)
  % Initialize latent variables
  F = [];
  covstruct = covstructfunc(theta_init);
  
  %get_rand_y = @get_rand_y_pcg;
  get_rand_y = @get_rand_y_gibbs;
  rand_y = get_rand_y(covstruct);

  func = @randfunc;
    function Y = randfunc(Y,Imv)
  
    if isempty(F)
      % Use some better initialization here..
      [Y, F] = rand_y(Y,Imv);
    end
    
    % whitened
    % V = do_something(covstruct, F);
    % [theta, covstruct] = rand_theta(Y, V);
    % F = do_it_back(covstruct, V);
    % conditioned
    % rand_theta(Y, F);
    
    % Sample hyperparameters
    [theta, covstruct_new] = rand_theta(Y, F);
    if ~isequal(covstruct, covstruct_new)
      covstruct = covstruct_new;
      rand_y = get_rand_y(covstruct);
      % If whitened likelihood, also update F..
    end
    
    % Sample Y
    [Y, F] = rand_y(Y,Imv);
    
    % Process samples??
    if ~isempty(samplefunc)
      samplefunc(Y, F, theta);
    end

    end

  end
  
  function [covstruct, covstruct_grad] = covstructfunc(theta, varargin)

  covstruct = [];
  covstruct_grad = [];
  
  theta1 = theta(1:n_theta1);
  theta2 = theta((n_theta1+1):(end-1));
  covstruct.s = theta(end);

  if nargout <= 1
    covstruct.K1 = covfunc1(theta1);
    covstruct.K2 = covfunc2(theta2);
  else
    [covstruct.K1, covstruct_grad.dK1] = covfunc1(theta1);
    [covstruct.K2, covstruct_grad.dK2] = covfunc2(theta2);
  end
  
  [linsolve1_f, logdet1_f] = get_solver1(covstruct.K1);
  [linsolve1_noisy, logdet1_noisy] = get_solver1(covstruct.K1, covstruct.s);

  [linsolve2_f, logdet2_f] = get_solver2(covstruct.K2);
  [linsolve2_noisy, logdet2_noisy] = get_solver2(covstruct.K2, covstruct.s);

  covstruct.linsolve_f = @(X) linsolve_kron(linsolve1_f, ...
                                            linsolve2_f, ...
                                            X);
  covstruct.logdet_f = @() (N2*logdet1_f()+N1*logdet2_f());
  
  covstruct.linsolve_z = @(X) linsolve_kron(linsolve1_noisy, ...
                                            linsolve2_noisy, ...
                                            X);
  covstruct.logdet_z = @() (N2*logdet1_noisy()+N1*logdet2_noisy());
  

  covstruct.multiply_f = @(X) kronprod(covstruct.K1,covstruct.K2,X);

  covstruct.multiply = @(X) (covstruct.multiply_f(X)+covstruct.s^2*X);
  multiply_pcg = @(X) reshape(covstruct.multiply(reshape(X, [N2 N1])), [N2*N1,1]);
  if true
    preconditioner_pcg = @(X) reshape(covstruct.linsolve_z(reshape(X, [N2 N1])), ...
                                      [N2*N1,1]);
  else
    preconditioner_pcg = [];
  end
  covstruct.linsolve = @(X) reshape(pcg(multiply_pcg, ...
                                        X(:), ...
                                        1e-4, ...
                                        100,  ...
                                        preconditioner_pcg, ...
                                        []), ...
                                    [N2, N1] );

  end
  
  function func = get_rand_y_pcg(covstruct)
  L1 = get_squareroot1(covstruct.K1);
  L2 = get_squareroot2(covstruct.K2);
  func = @rand_y_pcg;
    function [Y,F] = rand_y_pcg(Y, Imv)
    S = speye(numel(Y));
    S(Imv(:),:) = [];
    multiply = @(x) reshape(covstruct.multiply(reshape(x,size(Y))), ...
                            [numel(Y),1]);
    multiply_f = @(x) reshape(covstruct.multiply_f(reshape(x,size(Y))), ...
                              [numel(Y),1]);
    L1 = get_squareroot1(covstruct.K1);
    L2 = get_squareroot2(covstruct.K2);
    Z_f = gaussian_rand_kron(L1,L2);
    Z_noise = covstruct.s * randn(numel(Y),1);
    F = gaussian_rand_pcg(Y(:), ...
                          multiply, ...
                          multiply_f, ...
                          Z_f(:), ...
                          Z_noise(:), ...
                          [], ...
                          S, ...
                          'maxiter', 1e3, ...
                          'tol', 1e-5);
    F = reshape(F, size(Y));
    Y(Imv) = F(Imv) + covstruct.s*randn(size(Y(Imv)));
    end
  end

  function func = get_rand_y_gibbs(covstruct)

  L1 = get_squareroot1(covstruct.K1);
  L2 = get_squareroot2(covstruct.K2);
  func = @rand_y_gibbs;

    function [Y, F] = rand_y_gibbs(Y, Imv)
    % Sample latent function values
    F0 = gaussian_rand_kron(L1,L2);
    Y0 = F0 + covstruct.s*randn(N2,N1);
    F = F0 - covstruct.multiply_f( covstruct.linsolve( Y0 - Y ) );
    % Sample missing values
    Y(Imv) = F(Imv) + covstruct.s*randn(size(Y(Imv)));
    % You could evaluate:
    % rand(F), mean(F), mean(Y), 
    end
  
  end

  function func = get_loglikelihood_conditioned(covstruct, covstruct_grad)
  if nargin >= 2
    invK1 = sinv(covstruct.K1);
    invK2 = sinv(covstruct.K2);
  end
  func = @loglikelihood_conditioned;
    function [x, dx] = loglikelihood_conditioned(Y, F)
    Z_f = covstruct.linsolve_f(F);
    Y0 = (Y-F);
    x = gaussian_logpdf(F(:)'*Z_f(:), ...
                        0, ...
                        0, ...
                        covstruct.logdet_f(), ...
                        numel(F)) + ...
        gaussian_logpdf((Y0(:)'*Y0(:))/covstruct.s^2, ...
                        0, ...
                        0, ...
                        2*numel(Y)*log(covstruct.s), ...
                        numel(Y));
    if nargout >= 2
      n_theta = n_theta1 + n_theta2 + 1;
      dx = zeros(n_theta,1);
      % Gradient for K1 hyperparameters
      for n=1:n_theta1
        dx(n) = 0.5 * traceprod(Z_f,kronprod(covstruct_grad.dK1{n}, ...
                                             covstruct.K2, ...
                                             Z_f)) - ...
                0.5 * N2 * traceprod(invK1, ...
                                     covstruct_grad.dK1{n});
      end
      % Gradient for K2 hyperparameters
      for n=1:n_theta2
        n2 = n + n_theta1;
        dx(n2) = 0.5 * traceprod(Z_f,kronprod(covstruct.K1, ...
                                              covstruct_grad.dK2{n}, ...
                                              Z_f)) - ...
                 0.5 * N1 * traceprod(invK2, ...
                                      covstruct_grad.dK2{n});
      end
      % Gradient for noise parameter
      dx(end) = - 0.5 * traceprod(Y0,Y0) * (-2)*covstruct.s.^(-3) - ...
                numel(Y) / covstruct.s;
    end
    end
  end
  
  function func = get_loglikelihood_whitened_prior(covstruct, covstruct_grad)
  % marginalized
  % other options: whitened, conditioned
  func = @loglikelihood_whitened;
    function [x, dx] = loglikelihood_whitened_prior(Y, V)
    Z_v = covstruct.linsolve_z(V);
    error('not ready yet');
    %F =
    %Z_y = (Y-covstruct.);
    x = gaussian_logpdf(F(:)'*Z_f(:), ...
                        0, ...
                        0, ...
                        covstruct.logdet_f(), ...
                        numel(F)) + ...
        gaussian_logpdf((Z_y(:)'*Z_y(:))/covstruct.s^2, ...
                        0, ...
                        0, ...
                        2*numel(Y)*log(covstruct.s), ...
                        numel(Y));
    if nargout >= 2
      for n=1:N
        V = gradient{n}(Z);
        dx(n) = 0.5 * (V(:)'*Z(:)) - ...
                0.5 * trace_linsolve_gradient{n}();
        
      end
      error('hmm hmm');
    end
    end
  end
  
  function func = get_loglikelihood_whitened_what(covstruct, covstruct_grad)
  % marginalized
  % other options: whitened, conditioned
  func = @loglikelihood_whitened;
    function [x, dx] = loglikelihood_whitened(Y, V)
    Z_v = covstruct.linsolve_z(V);
    error('not ready yet');
    %F =
    %Z_y = (Y-covstruct.);
    x = gaussian_logpdf(F(:)'*Z_f(:), ...
                        0, ...
                        0, ...
                        covstruct.logdet_f(), ...
                        numel(F)) + ...
        gaussian_logpdf((Z_y(:)'*Z_y(:))/covstruct.s^2, ...
                        0, ...
                        0, ...
                        2*numel(Y)*log(covstruct.s), ...
                        numel(Y));
    if nargout >= 2
      for n=1:N
        V = gradient{n}(Z);
        dx(n) = 0.5 * (V(:)'*Z(:)) - ...
                0.5 * trace_linsolve_gradient{n}();
        
      end
      error('hmm hmm');
    end
    end
  end
  
end


function [linsolve_func, logdet_func] = get_covsolver_ldlchol(K, s2)
if nargin >= 2
  K = K + s2*speye(size(K));
end
LD = ldlchol(K);
linsolve_func = @(x) linsolve_ldlchol(LD,x);
logdet_func = @() logdet_ldlchol(LD);
end

function [linsolve_func, logdet_func] = get_covsolver_chol(K, s2)
if nargin >= 2
  K = K + s2*speye(size(K));
end
U = chol(K);
linsolve_func = @(x) linsolve_chol(U,x);
logdet_func = @() logdet_chol(U);
end

function [linsolve_func, logdet_func] = get_covsolver_lowrank(K, s2)
error('Not implemented yet');
end







































% $$$ function func = get_covstructfunc_kron3_pcg(covfunc1, covfunc2)
% $$$ 
% $$$ % additive noise:
% $$$ % cov = kron(K1, K2) + s^2
% $$$ 
% $$$ options = struct( ...
% $$$     'get_solver',         [],               ...
% $$$     'multiplier1',        @(K,x) (K*x),      ...
% $$$     'multiplier2',        @(K,x) (K*x),      ...
% $$$     'user_data',    []);
% $$$ 
% $$$ % Parse arguments
% $$$ [options, errmsg] = argparse( options, varargin{:} );
% $$$ error(errmsg);
% $$$ 
% $$$ % maybe this function should return the likelihood-function?
% $$$ 
% $$$ % solvers you might wanna consider:
% $$$ % - PCG (linsolve, logdet_noiseless, logdet_noise
% $$$ % - SVD/Schur (linsolve, logdet
% $$$ 
% $$$ Z1 = decompose(K1)
% $$$ Z2 = decompose(K2)
% $$$ @(x) linsolver(Z1, Z2, s, x)
% $$$ @(y) loglikelihood(Z1, Z2, s, y)
% $$$ 
% $$$ nin = nargin;
% $$$ if nargin == 2
% $$$   covfunc2 = get_solver1;
% $$$ end
% $$$ 
% $$$ [n_theta1, N1, M1] = covfunc1();
% $$$ [n_theta2, N2, M2] = covfunc2();
% $$$ func = @covstructfunc;
% $$$ 
% $$$   function covstruct = covstructfunc(theta)
% $$$   theta1 = theta(1:n_theta1);
% $$$   K1 = covfunc1(log(theta1));
% $$$   if nin > 2
% $$$     [linsolve1, logdet1] = get_solver1(K1);
% $$$     multiply1 = get_multiplier1(K1);
% $$$   else
% $$$     multiply1 = @(x) K1*x;
% $$$   end
% $$$ 
% $$$   theta2 = theta((n_theta1+1):(end-1));
% $$$   K2 = covfunc2(log(theta2));
% $$$   [linsolve2, logdet2] = get_solver2(K2);
% $$$   multiply2 = get_multiplier2(K2);
% $$$ 
% $$$   s = theta(end);
% $$$ 
% $$$   covstruct.multiply1 = multiply1;
% $$$   covstruct.multiply2 = multiply2;
% $$$   
% $$$   covstruct.multiply_noiseless = @(X) kronprod(multiply1,multiply2,X);
% $$$   %cov1 = @(x) (multiply1(x)+s1^2*x);
% $$$   %cov2 = @(x) (multiply2(x)+s2^2*x);
% $$$   covstruct.multiply = @(X) (kronprod(multiply1,multiply2,X) + s^2*X);
% $$$   preconditioner
% $$$   covstruct.linsolve = @(X) pcg(@(Y) covstruct.multiply(reshape(Y,size(X))), ...
% $$$                                 X(:), ...
% $$$                                 1e-3, ...
% $$$                                 100,  ...
% $$$                                 preconditioner);
% $$$                                 
% $$$   % covstruct.linsolve = @(X) linsolve_kron(linsolve1, linsolve2, X);
% $$$   covstruct.logdet = @() (N2*logdet1() + N1*logdet2());
% $$$   covstruct.N = N1*N2;
% $$$ 
% $$$   covstruct.get_rand = @get_rand;
% $$$   %covstruct.get_squareroot_K = @get_squareroot_K;
% $$$   %covstruct.get_squareroot_noise = @() (@(X) s1*s2*X);
% $$$     function [func_f, func_noise] = get_rand()
% $$$     L1 = get_squareroot1(K1);
% $$$     L2 = get_squareroot2(K2);
% $$$     func_f = @rand_f;
% $$$     func_noise = @() s1*s2*randn(N2,N1);
% $$$     if nargout < 2
% $$$       func_f = @() (func_f()+func_noise());
% $$$     end
% $$$       function f = rand_f()
% $$$       f = gaussian_rand_kron(L1,L2) + ...
% $$$           gaussian_rand_kron(L1,s2*speye(N2)) + ...
% $$$           gaussian_rand_kron(s1*speye(N1),L2);
% $$$       end
% $$$     end
% $$$     
% $$$   end
% $$$   
% $$$ 
% $$$ end

% $$$ function [x, dx] = loglikelihoodfunc(Y, covstruct)
% $$$ % marginalized
% $$$ % other options: whitened, conditioned
% $$$ Z = covstruct.linsolve(Y);
% $$$ x = gaussian_logpdf(...
% $$$     Y(:)'*Z(:), ...
% $$$     0, ...
% $$$     0, ...
% $$$     covstruct.logdet(), ...
% $$$     covstruct.N);
% $$$ if nargout >= 2
% $$$   for n=1:N
% $$$     V = covstruct.gradient{n}(Z);
% $$$     dx(n) = 0.5 * (V(:)'*Z(:)) - ...
% $$$             0.5 * covstruct.trace_linsolve_gradient{n}();
% $$$     
% $$$   end
% $$$   error('hmm hmm');
% $$$ end
% $$$ end

% $$$ function func = get_loglikelihoodfunc(covstruct)
% $$$ 
% $$$ logdet = covstruct.logdet();
% $$$ N = covstruct.N;
% $$$ func = @loglikelihoodfunc;
% $$$ 
% $$$   function loglikelihoodfunc(Y)
% $$$   Z = covstruct.linsolve(Y);
% $$$   logp = gaussian_logpdf(...
% $$$       Y(:)'*Z(:), ...
% $$$       0, ...
% $$$       0, ...
% $$$       logdet, ...
% $$$       N);
% $$$   end
% $$$ end

% $$$ function func = get_logposteriorfunc(logpriorfunc, loglikelihoodfunc)
% $$$ func = @logposteriorfunc;
% $$$   function logposterior = logposteriorfunc(Y, X, covstruct)
% $$$   logprior = logpriorfunc(X);
% $$$   loglikelihood = loglikelihoodfunc(Y, covstruct);
% $$$   logposterior = logprior + loglikelihood;
% $$$   %logpdf = funcsum(logprior, likelihood);
% $$$   end
% $$$ end

% $$$ function func = get_(logprior, get_loglikelihoodfunc)
% $$$ func = @get_logposteriorfunc;
% $$$   function logposterior = get_logposteriorfunc(X, Z)
% $$$   prior = logprior(X);
% $$$   likelihood = loglikelihood(Z);
% $$$   logpdf = funcsum(@(Y)prior, likelihood);
% $$$   end
% $$$ end













% $$$ function logp = gaussian_logpdf_kron(Y, linsolve_Cov1, linsolve_Cov2, ...
% $$$                                      logdet_Cov1, logdet_Cov2)
% $$$ % gaussian_logpdf_kron_ldlchol
% $$$ [N2,N1] = size(Y);
% $$$ logp = gaussian_logpdf(...
% $$$     Y(:)'*linsolve_kron(linsolve_Cov1, ...
% $$$                         linsolve_Cov2, ...
% $$$                         Y, ...
% $$$                         'vector'), ...
% $$$     0, ...
% $$$     0, ...
% $$$     N2*logdet_Cov1 + N1*logdet_Cov2, ...
% $$$     N1*N2);
% $$$ end

% $$$ function f(varargin)
% $$$   getfuncs = varargin;
% $$$   function g(varargin)
% $$$   funcs = cell(size(getfuncs));
% $$$   for n=1:numel(funcs)
% $$$     funcs{n} = getfuncs{n}(varargin{:});
% $$$   end
% $$$   h = @funcsum;
% $$$     function varargout = funcsum(varargin)
% $$$     for n=1:numel(funcs)
% $$$       outs;
% $$$     end
% $$$     end
% $$$   end
% $$$ end


% $$$ function f(f1, f2)
% $$$   function g(varargin)
% $$$   f1_eval = f1(varargin{:});
% $$$   f2_eval = f2(varargin{:});
% $$$   h = @(varargin) (f1_eval(varargin{:}) + f2_eval(varargin{:}));
% $$$   end
% $$$ end

% $$$ function what(who)
% $$$ theta_proposal = rand_proposal(theta);
% $$$ covstruct_proposal = covfunc(theta_proposal);
% $$$ loglikefunc_proposal = getfunc_loglikelihood(covstruct_proposal);
% $$$ loglike_proposal = loglikefunc_proposal(Y);
% $$$ logprior_proposal = logprior(theta);
% $$$ 
% $$$ ratio_proposal = (loglike_proposal + logprior_proposal) - ...
% $$$     logpdf_proposal(theta_proposal, theta);
% $$$ ratio_current = (loglike_current + logprior_current) - ...
% $$$     logpdf_proposal(theta, theta_proposal);
% $$$ r = ratio_proposal - ratio_current;
% $$$ 
% $$$ logproposal = logpdf_proposal(theta_proposal, theta) - ...
% $$$     logpdf_proposal(theta, theta_proposal);
% $$$ 
% $$$ loglike = loglikefunc_current(Y);
% $$$ 
% $$$ 
% $$$ end

% $$$ function p = mcmc_metropolishastings_acceptance(x_proposal, x, logdensity, ...
% $$$                                                 rand_proposal, logpdf_proposal)
% $$$ 
% $$$ r = (logdensity(x_proposal) - logpdf_proposal(x_proposal, x)) - ...
% $$$     (logdensity(x) - logpdf_proposal(x, x_proposal));
% $$$ p = (log(rand) < r)
% $$$ end
% $$$ 

% $$$ function func = get_logposterior(logprior, model)
% $$$ loglikelihood = model.get_loglikelihood();
% $$$ func = @logposterior;
% $$$   function l = logposterior(Y, X)
% $$$   l = logprior(X) + loglikelihood(Y);
% $$$   end
% $$$ end

% $$$ function func = mcmc_init_slicesampling(x_init, get_logpdf, func_x)
% $$$ % Axis aligned slice sampling
% $$$ x_current = x_init;
% $$$ fx_current = func_x(x_current);
% $$$ logpdf_current = get_logpdf(fx_current);
% $$$ x_proposal = x_current;
% $$$ func = @slicesampling;
% $$$   function [x, fx] = slicesampling(varargin)
% $$$   lp_current = logpdf_current(x_current, varargin{:});
% $$$   for i=1:numel(x_current)
% $$$     % Uniform vertical component
% $$$     lp_y = log(rand()) + lp_current;
% $$$     % Find the initial horizontal interval
% $$$     j = 0;
% $$$     while lp_min >= lp_current;
% $$$       xi_min = x_current(i) - 2^j;
% $$$       fx_min = func_x(xi_min, i, x_current, fx_current);
% $$$       logpdf_min = get_logpdf(fx_min);
% $$$       lp_min = logpdf_min(xi_min, varargin{:});
% $$$       j = j + 1;
% $$$     end
% $$$     j = 0;
% $$$     while lp_max >= lp_current;
% $$$       xi_max = x_current(i) + 2^j;
% $$$       fx_max = func_x(xi_max, i, x_current, fx_current);
% $$$       logpdf_max = get_logpdf(fx_max);
% $$$       lp_max = logpdf_max(xi_max, varargin{:});
% $$$       j = j + 1;
% $$$     end
% $$$     % Draw the horizontal component
% $$$     accepted = false;
% $$$     while ~accepted
% $$$       % Propose uniformly from the horizontal interval
% $$$       x_proposal(i) = unifrnd(xi_min, xi_max);
% $$$       % Get density value for the proposal
% $$$       fx_proposal = func_x(x_proposal, i, x_current, fx_current);
% $$$       logpdf_proposal = get_logpdf(fx_proposal);
% $$$       lp_proposal = logpdf_proposal(x_proposal, varargin{:});
% $$$       % Reject or accept the proposal
% $$$       if lp_proposal < lp_y
% $$$         % Reject the proposal and fix the horizontal interval
% $$$         if x_proposal(i) < x_current(i)
% $$$           xi_min = x_proposal(i);
% $$$         else
% $$$           xi_max = x_proposal(i);
% $$$         end
% $$$       else
% $$$         % Accept the proposal
% $$$         logpdf_current = logpdf_proposal;
% $$$         lp_current = lp_proposal;
% $$$         x_current = x_proposal;
% $$$         fx_current = fx_proposal;
% $$$         accepted = true;
% $$$       end
% $$$     end
% $$$   end
% $$$   x = x_current;
% $$$   fx = fx_current;
% $$$   end
% $$$ end

% $$$ function func = mcmc_init_metropolishastings(x_init, get_logpdf, q, logpdf_q, ...
% $$$                                              func_x)
% $$$ 
% $$$ if nargin < 6
% $$$   samplefunc = [];
% $$$ end
% $$$ x_current = x_init;
% $$$ fx_current = func_x(x_current);
% $$$ logpdf_current = get_logpdf(x_current, fx_current);
% $$$ 
% $$$ func = @metropolishastings;
% $$$ %model_init = mapx_current;
% $$$ 
% $$$ % $$$   [z, dz] = func(theta);
% $$$ % $$$   logpdf = get_logpdf(z, dz);
% $$$ % $$$   [lp, dlp] = logpdf(theta, Y);
% $$$ % $$$   [lp, dlp] = logpdf(theta, Y, V);
% $$$ % $$$   
% $$$   function [x, fx] = metropolishastings(varargin)
% $$$   x_proposal = q(x_current);
% $$$   fx_proposal = func_x(x_proposal);
% $$$   logpdf_proposal = get_logpdf(x_proposal, fx_proposal);
% $$$   lp_proposal = logpdf_proposal(varargin{:});
% $$$   % logpdf_proposal = logpdf(x_proposal, mapx_proposal);
% $$$   lp_current = logpdf_current(varargin{:});
% $$$   r = (lp_proposal - logpdf_q(x_proposal, x_current)) - ...
% $$$       (lp_current - logpdf_q(x_current, x_proposal));
% $$$   % r = (logpdf_proposal(Y) - logpdf_q(x_proposal, x)) - ...
% $$$   %     (logpdf_current(Y) - logpdf_q(x, x_proposal));
% $$$   if log(rand) < r
% $$$     disp('accept')
% $$$     fx_current = fx_proposal;
% $$$     x_current = x_proposal;
% $$$     logpdf_current = logpdf_proposal;
% $$$   else
% $$$     disp('reject')
% $$$   end
% $$$   x = x_current;
% $$$   fx = fx_current;
% $$$ % $$$   if ~isempty(samplefunc)
% $$$ % $$$     samplefunc(x);
% $$$ % $$$   end
% $$$   end
% $$$ end

% $$$ function [func, mapx_init] = mcmc_init_metropolishastings(x_init, logpdf, ...
% $$$                                                   q, logpdf_q, mapfunc, ...
% $$$                                                   samplefunc)
% $$$ 
% $$$ if nargin < 6
% $$$   samplefunc = [];
% $$$ end
% $$$ x_current = x_init;
% $$$ mapx_current = mapfunc(x_current);
% $$$ % logpdf_current = logpdf(x_current, mapx_current);
% $$$ 
% $$$ func = @metropolishastings;
% $$$ mapx_init = mapx_current;
% $$$ 
% $$$   function [x, mapx] = metropolishastings(Y)
% $$$   x_proposal = q(x_current);
% $$$   mapx_proposal = mapfunc(x_proposal);
% $$$   logpdf_proposal = logpdf(Y, x_proposal, mapx_proposal);
% $$$   % logpdf_proposal = logpdf(x_proposal, mapx_proposal);
% $$$   logpdf_current = logpdf(Y, x_current, mapx_current);
% $$$   r = (logpdf_proposal - logpdf_q(x_proposal, x_current)) - ...
% $$$       (logpdf_current - logpdf_q(x_current, x_proposal));
% $$$   % r = (logpdf_proposal(Y) - logpdf_q(x_proposal, x)) - ...
% $$$   %     (logpdf_current(Y) - logpdf_q(x, x_proposal));
% $$$   if log(rand) < r
% $$$     disp('accept')
% $$$     mapx_current = mapx_proposal;
% $$$     x_current = x_proposal;
% $$$     logpdf_current = logpdf_proposal;
% $$$   else
% $$$     disp('reject')
% $$$   end
% $$$   x = x_current;
% $$$   mapx = mapx_current;
% $$$   if ~isempty(samplefunc)
% $$$     samplefunc(x);
% $$$   end
% $$$   end
% $$$ end
% $$$ 

% $$$ function x = mcmc_metropolishastings_step(x, logdensity, rand_proposal, ...
% $$$                                           logpdf_proposal)
% $$$ 
% $$$ x_proposal = rand_proposal(x);
% $$$ r = (logdensity(x_proposal) - logpdf_proposal(x_proposal, x)) - ...
% $$$     (logdensity(x) - logpdf_proposal(x, x_proposal));
% $$$ if log(rand) < r
% $$$   x = x_proposal;
% $$$ end
% $$$ end

% $$$ function L = get_squareroot_chol(K)
% $$$ L = chol(K, 'lower');
% $$$ %squareroot_func = @(x) U'*x;
% $$$ end


%function get_covfunc(covfunc1,covfunc2,get_solver1,get_solver2)
%  function covmatrix = covfunc(theta)

% $$$ function covstruct = covstruct_kron_noisy(covfunc1, theta1, s1, get_solver1, ...
% $$$                                           covfunc2, theta2, s2, get_solver2)
% $$$ n_theta1 = covfunc1();
% $$$ theta1 = theta(1:n_theta1);
% $$$ K1 = covfunc1(log(theta1));
% $$$ Cov1 = K1 + s1*speye(size(K1));
% $$$ [linsolve1, logdet1] = get_solver1(Cov1);
% $$$ 
% $$$ n_theta2 = covfunc2();
% $$$ theta2 = theta((n_theta1+1):end);
% $$$ K2 = covfunc2(log(theta2));
% $$$ Cov2 = K1 + s2*speye(size(K2));
% $$$ [linsolve2, logdet2] = get_solver2(Cov2);
% $$$ 
% $$$ N1 = size(K1,1);
% $$$ N2 = size(K2,1);
% $$$ covstruct.K = @(X) (kronprod(K1,K2,X);
% $$$ covstruct.linsolve = @(X) linsolve_kron(linsolve1, linsolve2, X);
% $$$ covstruct.logdet = @() (N2*logdet1() + N1*logdet2());
% $$$ 
% $$$ end

% $$$ function covstruct = covstruct_kron(covfunc1, theta1, get_solver1, covfunc2, ...
% $$$                                     theta2, get_solver2)
% $$$ n_theta1 = covfunc1();
% $$$ theta1 = theta(1:n_theta1);
% $$$ K1 = covfunc1(log(theta1));
% $$$ [linsolve1, logdet1] = get_solver1(K1);
% $$$ 
% $$$ n_theta2 = covfunc2();
% $$$ theta2 = theta((n_theta1+1):end);
% $$$ K2 = covfunc2(log(theta2));
% $$$ [linsolve2, logdet2] = get_solver2(K2);
% $$$ 
% $$$ N1 = size(K1,1);
% $$$ N2 = size(K2,1);
% $$$ covstruct.K = @(X) kronprod(K1,K2,X);
% $$$ covstruct.linsolve = @(X) linsolve_kron(linsolve1, linsolve2, X);
% $$$ covstruct.logdet = @() (N2*logdet1() + N1*logdet2());
% $$$ 
% $$$ end

% $$$ function logpdffunc = get_gaussian_logpdf_kron(covmatrix)
% $$$ 
% $$$ logdet_Cov1 = covmatrix.logdet_Cov1();
% $$$ logdet_Cov2 = covmatrix.logdet_Cov2();
% $$$ logpdffunc = @(Y) gaussian_logpdf_kron(Y, ...
% $$$                                        covmatrix.linsolve_Cov1, ...
% $$$                                        covmatrix.linsolve_Cov2, ...
% $$$                                        logdet_Cov1, ...
% $$$                                        logdet_Cov2);
% $$$ end
% $$$ 
% $$$ function logp = gaussian_logpdf_kron(Y, linsolve_Cov1, linsolve_Cov2, ...
% $$$                                      logdet_Cov1, logdet_Cov2)
% $$$ % gaussian_logpdf_kron_ldlchol
% $$$ [N2,N1] = size(Y);
% $$$ logp = gaussian_logpdf(...
% $$$     Y(:)'*linsolve_kron(linsolve_Cov1, ...
% $$$                         linsolve_Cov2, ...
% $$$                         Y, ...
% $$$                         'vector'), ...
% $$$     0, ...
% $$$     0, ...
% $$$     N2*logdet_Cov1 + N1*logdet_Cov2, ...
% $$$     N1*N2);
% $$$ end


% $$$ function randfunc = get_randfunc_theta_metropolis(rand_proposal)
% $$$ 
% $$$   randfunc = @rand_theta;
% $$$ 
% $$$   function [covmatrix, theta] = rand_theta(theta, covmatrix, Y, covfunc, ...
% $$$                                            loglikelihood, logprior)
% $$$   % Sample new hyperparameters using Metropolis
% $$$   theta_prop = rand_proposal(theta);
% $$$   covmatrix_prop = covfunc(theta_prop);
% $$$   logposterior = loglikelihood(covmatrix,Y) + logprior(theta);
% $$$   logposterior_prop = loglikelihood(covmatrix_prop,Y) + logprior(theta_prop);
% $$$   if log(rand) < (logposterior_prop - logposterior)
% $$$     covmatrix = covmatrix_prop;
% $$$     theta = theta_prop;
% $$$     %  sample_variables = get_sampler();
% $$$     disp('accepted')
% $$$   else
% $$$     disp('rejected')
% $$$   end
% $$$   end
% $$$ end


% $$$ function [Y, F_noiseless, thetas] = rand_y(Y,I)
% $$$ F = sample_variables(Y, covmatrix);
% $$$ Y = sample_observations(Y, I, F
% $$$ 
% $$$ 
% $$$ 
% $$$ [F, F_noiseless] = sample_variables(Y);
% $$$ 
% $$$ % Reconstruct missing values
% $$$ if nargin < 2
% $$$   Y = F + covmatrix.s1*covmatrix.s2*randn(size(F));
% $$$ else
% $$$   Y(I) = F(I) + covmatrix.s1*covmatrix.s2*randn(size(F(I)));
% $$$ end
% $$$ thetas = theta;
% $$$ 
% $$$ Ymean = (n*Ymean + F) / (n+1);
% $$$ Ymean = (n*Ymean + F_noiseless) / (n+1);
% $$$ n = n+1;
% $$$ F_noiseless = Ymean;
% $$$ 
% $$$ %F_noiseless = F;
% $$$ 
% $$$ end
% $$$     function sample_hyperparameters(Y)
% $$$     % Sample new hyperparameters using Metropolis
% $$$     theta_prop = proprand(theta);
% $$$     covmatrix_prop = covfunc(theta_prop);
% $$$     if log(rand) < (targetdist(covmatrix_prop,Y) - targetdist(covmatrix,Y))
% $$$       covmatrix = covmatrix_prop;
% $$$       theta = theta_prop;
% $$$       sample_variables = get_sampler();
% $$$       disp('accepted')
% $$$     else
% $$$       disp('rejected')
% $$$     end
% $$$     end
% $$$ 
% $$$     function theta = proprand(theta)
% $$$     theta = exp(log(theta) + 0.001*randn(size(theta)));
% $$$     end
% $$$ 
% $$$   
% $$$     function logp = targetdist(covmat, Y)
% $$$     % gaussian_logpdf_kron_ldlchol
% $$$     logp = gaussian_logpdf(...
% $$$         Y(:)'*linsolve_kron_ldlchol(covmat.LD_Cov1,covmat.LD_Cov2,Y,'vector'), ...
% $$$         0, ...
% $$$         0, ...
% $$$         logdet_kron_ldlchol(covmat.LD_Cov1,covmat.LD_Cov2), ...
% $$$         numel(Y));
% $$$     end
% $$$ 
% $$$     function sampler_func = get_sampler()
% $$$     L1 = lchol(covmatrix.K1);
% $$$     L2 = lchol(covmatrix.K2);
% $$$     N1 = size(covmatrix.K1,1);
% $$$     N2 = size(covmatrix.K2,1);
% $$$     s1 = covmatrix.s1;
% $$$     s2 = covmatrix.s2;
% $$$     LD_Cov1 = covmatrix.LD_Cov1;
% $$$     LD_Cov2 = covmatrix.LD_Cov2;
% $$$     K1 = covmatrix.K1;
% $$$     K2 = covmatrix.K2;
% $$$     Covfun = @(Y) linsolve_kron_ldlchol(LD_Cov1,LD_Cov2,Y);
% $$$     Kfun_noisy = @(X) (kronprod(K1,K2,X) + ...
% $$$                        kronprod(K1,s2^2,X) + ...
% $$$                        kronprod(s1^2,K2,X));
% $$$     Kfun = @(X) kronprod(K1,K2,X);
% $$$     sampler_func = @sampler;
% $$$            
% $$$       function [F_noisy, F_noiseless] = sampler(Y)
% $$$       % function [F_noisy, F] = sampler(Y)
% $$$       z_k1k2 = kronprod(L1,L2,randn(N2,N1));
% $$$       z_k1s_sk2 = kronprod(L1,s2,randn(N2,N1)) + ...
% $$$           kronprod(s1,L2,randn(N2,N1));
% $$$       z_s = s1*s2*randn(N2,N1);
% $$$       F_noiseless = gaussian_rand_conditional(Y, Covfun, Kfun, z_k1k2, ...
% $$$                                               z_k1s_sk2+z_s);
% $$$       F_noisy = gaussian_rand_conditional(Y, Covfun, Kfun_noisy, z_k1k2+z_k1s_sk2, ...
% $$$                                           z_s);
% $$$       end
% $$$     end
% $$$ 
% $$$   end


% $$$ function func = get_logpdf_log(get_logpdf)
% $$$ func = @(logtheta, varargin) theta2logtheta(logtheta, ...
% $$$                                             get_logpdf(exp(logtheta), ...
% $$$                                                   varargin{:}));
% $$$   function func = theta2logtheta(logtheta, logpdf)
% $$$   func = @logpdf_transformed;
% $$$   z = exp(sum(logtheta));
% $$$     function [lp, dlp] = logpdf_transformed(varargin)
% $$$     lp = logpdf(varargin{:}) * z;
% $$$     end
% $$$   end
% $$$ end

% $$$ function init_maximumlikelihood(x_init, logpdf, covstructfunc)
% $$$ 
% $$$ x_current = x_init;
% $$$ 
% $$$   function [x, mapx] = ml(Y)
% $$$   [f, df_dx] = logpdf(Y, x, mapx_proposal)
% $$$   func = @(x) logpdf(Y, x, covstructfunc(x, 'gradient'))
% $$$   x_current = minimize(func, x_current);
% $$$   x = x_current;
% $$$   mapx = covstructfunc(x);
% $$$   end
% $$$   
% $$$ end

% $$$ function [covfunc, dcovfunc] = gp_cov_example(D)
% $$$   function K = cov(theta)
% $$$   K = theta(1) * exp(D.^2);
% $$$   end
% $$$   function dK = dcov(theta, K)
% $$$   dK = K / theta(1);
% $$$   end
% $$$ end
% $$$ 
% $$$ 
% $$$ function func = gp_cov_init_scale(cov)
% $$$ func = @gp_cov_scale;
% $$$   function [covfunc, dcovfunc] = gp_cov_scale(in)
% $$$     function cov(theta)
% $$$     end
% $$$     function dcov(theta, K)
% $$$     end
% $$$   end
% $$$ end

% $$$ function covstruct2 = gp_cov_scale2(covstruct1)
% $$$ covstruct2.cov = @cov;
% $$$ covstruct2.dcov = @dcov;
% $$$   function K = cov(theta)
% $$$   K = theta(1)^2 * covstruct1.cov(theta(2:end));
% $$$   % multiply = covstruct1.cov(theta(2:end));
% $$$   % multiply = @(x) (theta(1)^2 * multiply(x));
% $$$   end
% $$$   function dK = dcov(theta, K)
% $$$   dK = covstruct1.dcov(theta(2:end), K/theta(1)^2); %argh..
% $$$   dK = cat(3, 2/theta(1)*K, dK);
% $$$   end
% $$$ end
% $$$ 

% $$$ function samplefunc = get_sample_function(D_theta, N)
% $$$ results.Y = [];
% $$$ results.F1 = 0;
% $$$ results.F2 = 0;
% $$$ results.theta = zeros(D_theta, N);
% $$$ samplefunc = @process_sample;
% $$$   function res = process_sample(n, Y, theta, model)
% $$$   if nargin >= 1
% $$$     results.Y = Y;
% $$$     mu1 = model.cov.multiply_noiseless(model.cov.linsolve(Y));
% $$$     mu2 = kronprod(model.cov.multiply1, ...
% $$$                    model.cov.multiply2, ...
% $$$                    model.cov.linsolve(Y));
% $$$     switch 'sample'
% $$$ 
% $$$      case 'mean'
% $$$       results.F1 = (mu1 + (n-1)*results.F1) / n;
% $$$       results.F2 = (mu2 + (n-1)*results.F2) / n;
% $$$       
% $$$      case 'sample'
% $$$       results.F1 = mu1;
% $$$       results.F2 = mu2;
% $$$       
% $$$      case default
% $$$       error('Unknown storing');
% $$$       
% $$$     end
% $$$     results.theta(:,n) = theta(:);
% $$$   end
% $$$   if nargout >= 1
% $$$     res = results;
% $$$   end
% $$$   end
% $$$ end

% $$$ function gibbs(Y, I, N, rand_theta, samplefunc)
% $$$ % function gibbs(Y, I, covstruct, N)
% $$$ 
% $$$ % rand_y = get_rand_y(covstruct);
% $$$ 
% $$$ covstruct = [];
% $$$ 
% $$$ Y(I) = 0;
% $$$ 
% $$$ for n=1:N
% $$$   % Sample hyperparameters
% $$$   [theta, covstruct_new] = rand_theta(Y);
% $$$   if ~isequal(covstruct, covstruct_new)
% $$$     covstruct = covstruct_new;
% $$$     rand_y = covstruct.get_rand_y();
% $$$     %rand_y = get_rand_y(covstruct);
% $$$   end
% $$$ 
% $$$   % Sample Y
% $$$   Y = rand_y(Y,I);
% $$$   
% $$$   % Process samples
% $$$   samplefunc(n, Y, theta, covstruct);
% $$$   
% $$$ end
% $$$ 
% $$$ end


% $$$ function func = get_rand_y_gibbs(covstruct)
% $$$ 
% $$$ [rand_f,rand_noise] = covstruct.get_rand();
% $$$ func = @rand_y_gibbs;
% $$$ 
% $$$   function Y = rand_y_gibbs(Y, Imv)
% $$$   % Sample latent function values
% $$$   z_f = rand_f();
% $$$   z_y = z_f + rand_noise();
% $$$   F = z_f - covstruct.multiply_noiseless( covstruct.linsolve( z_y - Y ) );
% $$$   % Sample missing values
% $$$   z_noise = rand_noise();
% $$$   Y(Imv) = F(Imv) + z_noise(Imv);
% $$$   end
% $$$   
% $$$ end
% $$$ 
% $$$ function func = get_rand_y_pcg(covstruct)
% $$$ 
% $$$ rand_y = covstruct.get_rand();
% $$$ func = @rand_y_pcg;
% $$$ 
% $$$   function Y = rand_y_pcg(Y, Imv)
% $$$   N = numel(Y);
% $$$   I = speye(N);
% $$$   S = I(~Imv(:),:);
% $$$   %cov1 = @(x) covstruct.multiply1(
% $$$   multiply = @(x) reshape(covstruct.multiply(reshape(x,size(Y))),[N,1]);
% $$$   Z = rand_y();
% $$$   Yh = gaussian_rand_pcg(Y(:), ...
% $$$                          multiply, ...
% $$$                          multiply, ...
% $$$                          Z(:), ...
% $$$                          0, ...
% $$$                          [], ... % covstruct.linsolve
% $$$                          S, ...
% $$$                          'tol', 1e-3);
% $$$   Y(Imv) = Yh(Imv);
% $$$   end
% $$$   
% $$$ end


% $$$ function func = get_covstructfunc_kron2(covfunc1, get_solver1, get_multiplier1, ...
% $$$                                                   get_squareroot1, covfunc2, ...
% $$$                                                   get_solver2, get_multiplier2, ...
% $$$                                                   get_squareroot2)
% $$$ % additive noise to both kron-components:
% $$$ % kron(K1+s1^2, K2+s2^2)
% $$$ 
% $$$ nin = nargin;
% $$$ if nargin == 2
% $$$   covfunc2 = get_solver1;
% $$$ end
% $$$ 
% $$$ [n_theta1, N1, M1] = covfunc1();
% $$$ [n_theta2, N2, M2] = covfunc2();
% $$$ func = @covstructfunc;
% $$$ 
% $$$   function model = covstructfunc(theta)
% $$$   theta1 = theta(1:n_theta1);
% $$$   s1 = theta(n_theta1+1);
% $$$   K1 = covfunc1(theta1);
% $$$   if nin > 2
% $$$     [linsolve1, logdet1] = get_solver1(K1, s1^2);
% $$$     multiply1 = get_multiplier1(K1);
% $$$   else
% $$$     multiply1 = @(x) K1*x;
% $$$   end
% $$$ 
% $$$   theta2 = theta((n_theta1+2):(end-1));
% $$$   s2 = theta(end);
% $$$   K2 = covfunc2(theta2);
% $$$   [linsolve2, logdet2] = get_solver2(K2, s2^2);
% $$$   multiply2 = get_multiplier2(K2);
% $$$ 
% $$$   cov.multiply1 = multiply1;
% $$$   cov.multiply2 = multiply2;
% $$$   
% $$$   cov.multiply_noiseless = @(X) ( kronprod(multiply1,multiply2,X) + ...
% $$$                                   kronprod(multiply1,s2^2,X) + ...
% $$$                                   kronprod(s1^2,multiply2,X) );
% $$$   cov1 = @(x) (multiply1(x)+s1^2*x);
% $$$   cov2 = @(x) (multiply2(x)+s2^2*x);
% $$$   cov.multiply = @(X) (kronprod(cov1,cov2,X));
% $$$   cov.linsolve = @(X) linsolve_kron(linsolve1, linsolve2, X);
% $$$   cov.logdet = @() (N2*logdet1() + N1*logdet2());
% $$$   cov.N = N1*N2;
% $$$ 
% $$$   model.get_rand_y = @get_rand_y_gibbs;
% $$$   model.get_loglikelihood = @get_loglikelihood_marginalized;
% $$$   model.cov = cov;
% $$$ 
% $$$     function func = get_rand_y_gibbs()
% $$$ 
% $$$     %[rand_f,rand_noise] = covstruct.get_rand();
% $$$     L1 = get_squareroot1(K1);
% $$$     L2 = get_squareroot2(K2);
% $$$     func = @rand_y_gibbs;
% $$$ 
% $$$       function Y = rand_y_gibbs(Y, Imv)
% $$$       % Sample latent function values
% $$$       F0 = gaussian_rand_kron(L1,L2) + ...
% $$$            gaussian_rand_kron(L1,s2*speye(N2)) + ...
% $$$            gaussian_rand_kron(s1*speye(N1),L2);
% $$$       Y0 = F0 + s1*s2*randn(N2,N1);
% $$$       F = F0 - cov.multiply_noiseless( cov.linsolve( Y0 - Y ) );
% $$$       % Sample missing values
% $$$       Y(Imv) = F(Imv) + s1*s2*randn(size(Y(Imv)));
% $$$       end
% $$$   
% $$$     end
% $$$ 
% $$$     function func = get_loglikelihood_marginalized()
% $$$     % marginalized
% $$$     % other options: whitened, conditioned
% $$$     ld = cov.logdet();
% $$$     func = @loglikelihood_marginalized;
% $$$       function [x, dx] = loglikelihood_marginalized(Y)
% $$$       Z = cov.linsolve(Y);
% $$$       x = gaussian_logpdf(...
% $$$           Y(:)'*Z(:), ...
% $$$           0, ...
% $$$           0, ...
% $$$           ld, ...
% $$$           cov.N);
% $$$       if nargout >= 2
% $$$         for n=1:N
% $$$           V = gradient{n}(Z);
% $$$           dx(n) = 0.5 * (V(:)'*Z(:)) - ...
% $$$                   0.5 * trace_linsolve_gradient{n}();
% $$$           
% $$$         end
% $$$         error('hmm hmm');
% $$$       end
% $$$       end
% $$$     end
% $$$ 
% $$$   end
% $$$   
% $$$ end


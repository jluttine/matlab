function [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron1(covfunc1, solver1, covfunc2, solver2, logprior, dlogprior, ...
                  samplefunc, varargin)
% $$$ function [func1, func2] = gp_init_kron1(covfunc1, ...
% $$$                                         get_solver1,     ...
% $$$                                         get_squareroot1, ...
% $$$                                         covfunc2,        ...
% $$$                                         get_solver2,     ...
% $$$                                         get_squareroot2, ...
% $$$                                         logprior,        ...
% $$$                                         samplefunc)
% additive noise to both kron-components:
% kron(K1+s1^2, K2+s2^2)

[n_theta1, N1] = covfunc1();
[n_theta2, N2] = covfunc2();

options = struct( ...
    'noise_scale1', ones(N1,1), ...
    'noise_scale2', ones(N2,1));

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Sampler for theta
get_loglikelihood = @get_loglikelihood_marginalized;
get_dloglikelihood = @get_dloglikelihood_marginalized;

get_logpdf = @get_logposterior;
get_dlogpdf = @get_dlogposterior;
get_rand_model = @get_randfunc;
func_theta = @covstructfunc;

  function logpost = get_logposterior(f_theta)
  lprior = logprior(f_theta.theta);
  loglike = get_loglikelihood(f_theta);
  logpost = @logposterior;
    function lpost = logposterior(varargin)
    llike = loglike(varargin{:});
    lpost = lprior + llike;
    end
  end

  function dlogpost = get_dlogposterior(df_theta)
  dlprior = dlogprior(df_theta.theta);
  dloglike = get_dloglikelihood(df_theta);
  dlogpost = @dlogposterior;
    function dlpost = dlogposterior(varargin)
    dllike = dloglike(varargin{:});
    dlpost = dlprior + dllike;
    end
  end

  function func = get_randfunc(rand_theta, covstruct_init)
  % Initialize latent variables
  F = [];%zeros(N2,N1);
  covstruct = covstruct_init;
  %covstruct = covstructfunc(theta_init);
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
    [~, covstruct_new] = rand_theta(Y);
    if ~isequal(covstruct, covstruct_new)
      covstruct = covstruct_new;
      rand_y = get_rand_y(covstruct);
    end
    
    % Sample Y
    [Y, F] = rand_y(Y,I);
    
    % Process samples??
    if ~isempty(samplefunc)
      samplefunc(Y, F, covstruct.theta);
    end
    end
  end

  function [covstruct, dcovstruct] = covstructfunc(theta, varargin)
  
  % Parameters: a, theta1, s1, theta2, s2
  
  % a * (k_1 + s_1) * (k_2 + s_2)
  
  ind1 = 1 + (1:n_theta1);
  ind2 = (1 + n_theta1 + 1) + (1:n_theta2);
  
  theta1 = theta(ind1);
  theta2 = theta(ind2);

  if nargin == 1

    covstruct.a = theta(1);
    covstruct.s1 = theta(1+n_theta1+1);
    covstruct.s2 = theta(1+n_theta1+1+n_theta2+1);
    covstruct.S1 = spdiag(covstruct.s1*options.noise_scale1);
    covstruct.S2 = spdiag(covstruct.s2*options.noise_scale2);
    
    covstruct.theta = theta;
    
    if nargout <= 1
      K1 = covfunc1(theta1);
      K2 = covfunc2(theta2);
    else
      [K1, dK1] = covfunc1(theta1);
      [K2, dK2] = covfunc2(theta2);
    end

    covstruct.K1 = K1;
    covstruct.K2 = K2;

    %size(covstruct.K2)
    %size(covstruct.S2)
    covstruct.L1 = solver1.decompose(covstruct.K1+covstruct.S1.^2);
    covstruct.L2 = solver2.decompose(covstruct.K2+covstruct.S2.^2);

    if nargout >= 2
      dcovstruct = covstruct;
      dcovstruct.dK1 = dK1;
      dcovstruct.dK2 = dK2;
    end

  else
    
    if nargout > 1
      error('Not implemented for getting the gradient');
    end
    
    ind = varargin{1};
    covstruct_old = varargin{2};
    
    covstruct = covstruct_old;
    
    covstruct.a = theta(1);
    covstruct.s1 = theta(1+n_theta1+1);
    covstruct.s2 = theta(1+n_theta1+1+n_theta2+1);
    covstruct.theta = theta;
    covstruct.S1 = spdiag(covstruct.s1*options.noise_scale1);
    covstruct.S2 = spdiag(covstruct.s2*options.noise_scale2);
    
    if any(ind > 1)
      if any(ind <= (n_theta1+2))
        %disp('Update 1')
        covstruct.K1 = covfunc1(theta1);
        covstruct.L1 = solver1.decompose(covstruct.K1 + covstruct.S1.^2);
        %covstruct.L1 = solver1.decompose(covstruct.K1+covstruct.s1^2*speye(N1));
      end
      if any(ind >= (n_theta1+2))
        %disp('Update 2')
        covstruct.K2 = covfunc2(theta2);
        covstruct.L2 = solver2.decompose(covstruct.K2 + covstruct.S2.^2);
% $$$         covstruct.L2 = solver2.decompose(covstruct.K2 + ...
% $$$                                          covstruct.s2^2*speye(N2));
      end
    end
    
  end
  
  end

  function func = get_rand_y_gibbs(covstruct)

  multiply_noiseless = @(X) (kronprod(covstruct.K1,covstruct.K2,X) + ...
                             kronprod(covstruct.K1,covstruct.S2.^2,X) + ...
                             kronprod(covstruct.S1.^2,covstruct.K2,X));
  multiply_f = @(X) (kronprod(covstruct.K1,covstruct.K2,X));
  linsolve = @(X) linsolve_kron(@(x1) solver1.linsolve(covstruct.L1, x1), ...
                                @(x2) solver2.linsolve(covstruct.L2, x2), ...
                                X);

  LD = solver1.decompose(covstruct.K1);
  L1 = solver1.squareroot(LD);
  LD = solver2.decompose(covstruct.K2);
  L2 = solver2.squareroot(LD);
% $$$   [L1,p] = chol(covstruct.K1, 'lower');
% $$$   if p>0
% $$$     theta = covstruct.theta
% $$$     error('Matrix mus be positive definite.');
% $$$   end
% $$$   [L2,p] = chol(covstruct.K2, 'lower');
% $$$   if p>0
% $$$     theta = covstruct.theta
% $$$     error('Matrix mus be positive definite.');
% $$$   end

  func = @rand_y_gibbs;

    function [Y, F] = rand_y_gibbs(Y, Imv)
    % Sample latent function values
    F0 = covstruct.a * gaussian_rand_kron(L1,L2);
    F0_1 = covstruct.a * gaussian_rand_kron(L1,covstruct.S2);
    F0_2 = covstruct.a * gaussian_rand_kron(covstruct.S1,L2);
    Y0 = F0 + F0_1 + F0_2 + covstruct.a*kronprod(covstruct.S1,covstruct.S2,randn(N2,N1));
    Z = linsolve( Y0 - Y );
    F = F0 + F0_1 + F0_2 - multiply_noiseless(Z);
    % Sample missing values
    Z = covstruct.a*kronprod(covstruct.S1,covstruct.S2,randn(N2,N1));
    Y(Imv) = F(Imv) + Z(Imv);
    %Y(Imv) = F(Imv) + covstruct.a*covstruct.s1*covstruct.s2*randn(size(Y(Imv)));
    % Interesting function values
    %F = F0 - multiply_f(Z);
    % You could evaluate:
    % rand(F), mean(F), mean(Y), 
    end
  
  end

  function func = get_rand_y_pcg(covstruct)
  multiply_f = @(X) (kronprod(covstruct.K1,covstruct.K2,X));
  multiply_noiseless = @(X) (kronprod(covstruct.K1,covstruct.K2,X) + ...
                             kronprod(covstruct.K1,covstruct.S2.^2,X) + ...
                             kronprod(covstruct.S1.^2,covstruct.K2,X));
  multiply = @(X) (multiply_noiseless(X) + ...
                   kronprod(covstruct.S1.^2, covstruct.S2.^2, X));
  LD = solver1.decompose(covstruct.K1);
  L1 = solver1.squareroot(LD);
  LD = solver2.decompose(covstruct.K2);
  L2 = solver2.squareroot(LD);
% $$$   L1 = chol(covstruct.K1, 'lower');
% $$$   L2 = chol(covstruct.K2, 'lower');
  func = @rand_y_pcg;

    function [Y,F] = rand_y_pcg(Y, Imv)
% $$$     S = speye(numel(Y));
% $$$     S(Imv(:),:) = [];
% $$$     multiply = @(x) reshape(multiply(reshape(x,size(Y))), ...
% $$$                             [numel(Y),1]);
% $$$     multiply_noiseless = @(x) reshape(multiply_noiseless(reshape(x,size(Y))), ...
% $$$                               [numel(Y),1]);
    F = covstruct.a * gaussian_rand_kron(L1,L2) ...
        + covstruct.a * gaussian_rand_kron(L1,covstruct.S2) ...
        + covstruct.a * gaussian_rand_kron(covstruct.S1,L2);
    %Y0 = F0 + F0_1 + F0_2 + covstruct.s1*covstruct.s2*randn(N2,N1);
    %Z_f = gaussian_rand_kron(L1,L2);
    Z_noise = covstruct.a * kronprod(covstruct.S1, covstruct.S2, randn(size(Y)));
    %Z_noise = covstruct.a*covstruct.s1*covstruct.s2 * randn(numel(Y),1);
    % TODO: IS BELOW CORRECT WRT COVSTRUCT.A ???
    F = gaussian_rand_conjgradmv(multiply_noiseless, ...
                                 multiply, ...
                                 F, ...
                                 F + Z_noise - Y, ...
                                 ~Imv, ...
                                 'maxiter', 1e3, ...
                                 'tol', 1e-6, ...
                                 'verbose', true);
% $$$     F = gaussian_rand_pcg(Y(:), ...
% $$$                           multiply, ...
% $$$                           multiply_noiseless, ...
% $$$                           F0(:) + F0_1(:) + F0_2(:), ...
% $$$                           Z_noise(:), ...
% $$$                           [], ...
% $$$                           S, ...
% $$$                           'maxiter', 1e3, ...
% $$$                           'tol', 1e-3);
%    F = reshape(F, size(Y));
    Z_noise = covstruct.a*kronprod(covstruct.S1,covstruct.S2,randn(N2,N1));
    Y(Imv) = F(Imv) + Z_noise(Imv);
    %Y(Imv) = F(Imv) + covstruct.a*covstruct.s1*covstruct.s2*randn(size(Y(Imv)));
% $$$     multiply_f = @(x) reshape(multiply_f(reshape(x,size(Y))), ...
% $$$                               [numel(Y),1]);
% $$$     F = gaussian_rand_pcg(Y(:), ...
% $$$                           multiply, ...
% $$$                           multiply_f, ...
% $$$                           F0(:) + F0_1(:) + F0_2(:), ...
% $$$                           Z_noise(:), ...
% $$$                           [], ...
% $$$                           S, ...
% $$$                           'maxiter', 1e3, ...
% $$$                           'tol', 1e-3);
    end
  end

  function func = get_loglikelihood_marginalized(covstruct)
  % marginalized
  linsolve = @(X) linsolve_kron(@(x1) solver1.linsolve(covstruct.L1, x1), ...
                                @(x2) solver2.linsolve(covstruct.L2, x2), ...
                                X);
  ldet = N2*solver1.logdet(covstruct.L1) + ...
         N1*solver2.logdet(covstruct.L2) + ...
         2*N2*N1*log(covstruct.a);
  func = @loglikelihood_marginalized;
    function x = loglikelihood_marginalized(Y)
    Z = linsolve(Y);
    x = gaussian_logpdf(...
        (Y(:)'*Z(:))/covstruct.a^2, ...
        0, ...
        0, ...
        ldet, ...
        N1*N2);
    end
  end
  function func = get_dloglikelihood_marginalized(dcovstruct)
  % marginalized
  invCov1 = solver1.inv(dcovstruct.L1);
  invCov2 = solver2.inv(dcovstruct.L2);
  linsolve = @(X) linsolve_kron(@(x1) solver1.linsolve(dcovstruct.L1, x1), ...
                                @(x2) solver2.linsolve(dcovstruct.L2, x2), ...
                                X);
  %ldet = N2*solver1.logdet(dcovstruct.L1) + N1*solver2.logdet(dcovstruct.L2);
  func = @dloglikelihood_marginalized;
    function dx = dloglikelihood_marginalized(Y)
    Z = linsolve(Y);

    n_theta = n_theta1 + n_theta2 + 3;
    dx = zeros(n_theta,1);
    % Gradient for overall magnitude
    dx(1) = - 0.5 * traceprod(Y, Z) * (-2)*dcovstruct.a^(-3) ...
            - N1*N2 / dcovstruct.a;
    %dx(1) = 0.5 * traceprod(Z,kronprod(dK, Cov2, Z)) - ...
            %0.5 * N2 * traceprod(invCov1, dK);
    % Gradient for Cov1 hyperparameters
    Cov2 = dcovstruct.K2 + dcovstruct.S2.^2;
    for n=1:n_theta1
      n1 = 1 + n;
      dx(n1) = 0.5 * traceprod(Z,kronprod(dcovstruct.dK1{n}, ...
                                          Cov2, ...
                                          Z)) / dcovstruct.a^2 - ...
               0.5 * N2 * traceprod(invCov1, ...
                                    dcovstruct.dK1{n});
    end
    dx(1+n_theta1+1) = 0.5 * traceprod(Z,kronprod(2*dcovstruct.S1.^2/dcovstruct.s1, ...
                                                  Cov2, ...
                                                  Z)) / dcovstruct.a^2 - ...
        0.5 * N2 * sum(diag(invCov1).*diag(dcovstruct.S1).^2) * 2 / dcovstruct.s1;
    % Gradient for Cov2 hyperparameters
    Cov1 = dcovstruct.K1 + dcovstruct.S1.^2;
% $$$     Cov1 = dcovstruct.K1 + dcovstruct.s1^2*speye(N1);
    for n=1:n_theta2
      n2 = 1 + n_theta1 + 1 + n;
      dx(n2) = 0.5 * traceprod(Z,kronprod(Cov1, ...
                                          dcovstruct.dK2{n}, ...
                                          Z)) / dcovstruct.a^2- ...
               0.5 * N1 * traceprod(invCov2, ...
                                    dcovstruct.dK2{n});
    end
    dx(end) = 0.5 * traceprod(Z,kronprod(Cov1, ...
                                         2*dcovstruct.S2.^2/dcovstruct.s2, ...
                                         Z)) / dcovstruct.a^2 - ...
              0.5 * N1 * sum(diag(invCov2).*diag(dcovstruct.S2).^2) * 2 / dcovstruct.s2;
% $$$               0.5 * N1 * sum(diag(invCov2)) * 2*dcovstruct.s2;
    end
  end
  
end



% K = SUM_d a_d^2*(K{1,d} * K{2,d}) + s^2*DIAG(W.^2)
%
% THETA = [A_1, THETA_11, THETA_21, A_2, THETA_12, THETA_22, ..., S]
%
% Make it such that DIAG(W) must be a Kronecker product of two diagonal
% matrices DIAG(W1) and DIAG(W2). These are fixed "spatial" and "temporal"
% noise scales.

% This version can subtract approximate posterior mean from F when doing
% whitening for hyperparameter sampling.

function [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron6(covfuncs, solvers, logprior, dlogprior, varargin)


options = struct( ...
    'samplefunc',  [], ...
    'likelihood',  'whitened_prior', ...
    'noise_scale1', [], ...
    'noise_scale2', [], ...
    'rand_y',      'pcg'); % gibbs / pcg

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Number of Kronecker-product covariance functions in the sum
[D_kron, D_sum] = size(covfuncs);
if D_kron ~= 2
  error(['Supporting only products over two domains, that is, COVFUNCS ' ...
         'must have exactly two rows.'])
end

% Check the dimensionalities of the parameter vector for each covariance
% function and compute the indices for each in the joint parameter vector.
N = nan(D_kron,1);
M = nan(D_kron,1);
n_theta = nan(D_kron, D_sum);
ind_theta = cell(D_kron, D_sum);
ind_a = nan(D_sum,1);
ind_s = NaN;
iter = 1;
for j=1:D_sum
  ind_a(j) = iter;
  iter = iter + 1;
  for i=1:D_kron
    %[n_theta(i,j),~,~] = feval(covfuncs{i,j});
    [n_theta(i,j),N(i),M(i)] = feval(covfuncs{i,j});
    ind_theta{i,j} = iter + (1:n_theta(i,j)) - 1;
    iter = iter + n_theta(i,j);
  end
end
ind_s = iter;

% 2-D Kronecker hard coded here..
w = cell(D_kron,1);
if isempty(options.noise_scale1)
  w{1} = ones(N(1),1);
else
  w{1} = reshape(options.noise_scale1, N(1), 1);
end
if isempty(options.noise_scale2)
  w{2} = ones(N(2),1);
else
  w{2} = reshape(options.noise_scale2, N(2), 1);
end
W = w{2}*w{1}';

debug = false;

% Posterior density function
switch options.likelihood
% $$$  case 'conditioned'
% $$$   get_loglikelihood = @get_loglikelihood_conditioned;
% $$$   get_dloglikelihood = @get_dloglikelihood_conditioned;
 case 'whitened_prior'
  get_loglikelihood = @get_loglikelihood_whitened_prior;
  get_dloglikelihood = @get_dloglikelihood_whitened_prior;
 case 'whitened_posterior'
  get_loglikelihood = @get_loglikelihood_whitened_posterior;
  get_dloglikelihood = @get_dloglikelihood_whitened_posterior;
 case 'whitened_posterior2'
  get_loglikelihood = @get_loglikelihood_whitened_posterior2;
  get_dloglikelihood = @get_dloglikelihood_whitened_posterior2;
 case 'whitened_posterior3'
  get_loglikelihood = @get_loglikelihood_whitened_posterior3;
  get_dloglikelihood = @get_dloglikelihood_whitened_posterior3;
 otherwise
  error('Unknown likelihood')
end

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
  F = [];
  
  s2W2 = NaN;

  covstruct = covstruct_init;
  
  switch options.rand_y
   case 'gibbs'
    get_rand_y = @get_rand_y_gibbs;
   case 'pcg'
    get_rand_y = @get_rand_y_pcg;
  end
  rand_y = get_rand_y(covstruct);

  func = @randfunc;
    function Y = randfunc(Y,Imv)
  
    %    if isempty(F)
    % Use some better initialization here..
    %      fprintf('Draw Y and F..\n');
    t = cputime();
    [Y, F] = rand_y(Y,Imv);
      %fprintf(' %f seconds\n', cputime()-t);
%    end
    
%    fprintf('Draw theta..\n');
%    t = cputime();

    switch options.likelihood
      
     case 'whitened_prior'
      % Transform F to V which has a whitened prior
% $$$       U = cell(D_kron, D_sum);
% $$$       linsolver = cell(D_kron, D_sum);
      Z = zeros(size(F));
      for j=1:D_sum
        % Compute chol(K)'\F by chol(K)*(K\F) in order to utilize the
        % linsolve-functions
        Z(:,:,j) = linsolve_kron(covstruct.L{1,j}, ...
                                 covstruct.L{2,j}, ...
                                 F(:,:,j)) / covstruct.a(j);
      end
      
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, Z, Imv);
      
      if ~isequal(covstruct, covstruct_new)
        covstruct = covstruct_new;
        rand_y = get_rand_y(covstruct);
      end
      
      % Transform whitened variable V to the model variable F
      for j=1:D_sum
        F(:,:,j) = covstruct.a(j) * kronprod(covstruct.L{1,j}, ...
                                             covstruct.L{2,j}, ...
                                             Z(:,:,j));
      end
      
     case 'whitened_posterior'

      % Transform F to Z
      a = [covstruct.a(:); covstruct.s];
      V = linsolve_kron(@(x) feval(solvers{1,1}.linsolve, ...
                                   covstruct.LD_Cov{1}, ...
                                   x), ...
                        @(x) feval(solvers{2,1}.linsolve, ...
                                   covstruct.LD_Cov{2}, ...
                                   x), ...
                        Y) * sum(a)^2 / sum(a.^2);
      Z = zeros(size(F));
      for j=1:D_sum
        % Approximate posterior mean
        Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                         covstruct.K{2,j}, ...
                                         V);
        % Approximate whitening
        Z(:,:,j) = linsolve_kron(covstruct.L{1,j}, ...
                                 covstruct.L{2,j}, ...
                                 F(:,:,j) - Fh) / covstruct.a(j);
      end
      
      %error('debugging')
    
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, Z, Imv);
    
      if ~isequal(covstruct, covstruct_new)
        covstruct = covstruct_new;
        rand_y = get_rand_y(covstruct);
      end
      
      % Transform Z back to F
      a = [covstruct.a(:); covstruct.s];
      V = linsolve_kron(@(x) feval(solvers{1,1}.linsolve, ...
                                   covstruct.LD_Cov{1}, ...
                                   x), ...
                        @(x) feval(solvers{2,1}.linsolve, ...
                                   covstruct.LD_Cov{2}, ...
                                   x), ...
                        Y) * sum(a)^2 / sum(a.^2);
      for j=1:D_sum
        % Approximate posterior mean (bias term)
        Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                         covstruct.K{2,j}, ...
                                         V);
        % Unwhiten Z
        X = kronprod(covstruct.L{1,j}, ...
                     covstruct.L{2,j}, ...
                     Z(:,:,j)) * covstruct.a(j);
        % Sum the unwhitened Z and the bias
        F(:,:,j) = Fh + X;
      end
      
     case 'whitened_posterior2'
      
      % Transform F to Z
      s2W2 = covstruct.s^2*W.^2;
      V = conjgradmv(@multiply_y, ...
                     Y, ...
                     ~Imv, ...
                     'maxiter', 1e3, ...
                     'tol', 1e-6, ...
                     'verbose', true);
      Z = zeros(size(F));
      for j=1:D_sum
        % Approximate posterior mean
        Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                         covstruct.K{2,j}, ...
                                         V);
        % Approximate whitening
        Z(:,:,j) = linsolve_kron(covstruct.L{1,j}, ...
                                 covstruct.L{2,j}, ...
                                 F(:,:,j) - Fh) / covstruct.a(j);
      end
      
      %error('debugging')
    
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, Z, Imv);
    
      if ~isequal(covstruct, covstruct_new)
        covstruct = covstruct_new;
        rand_y = get_rand_y(covstruct);
      end
      
      % Transform Z back to F
      s2W2 = covstruct.s^2*W.^2;
      V = conjgradmv(@multiply_y, ...
                     Y, ...
                     ~Imv, ...
                     'maxiter', 1e3, ...
                     'tol', 1e-6, ...
                     'verbose', true);
      for j=1:D_sum
        % Approximate posterior mean (bias term)
        Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                         covstruct.K{2,j}, ...
                                         V);
        % Unwhiten Z
        X = kronprod(covstruct.L{1,j}, ...
                     covstruct.L{2,j}, ...
                     Z(:,:,j)) * covstruct.a(j);
        % Sum the unwhitened Z and the bias
        F(:,:,j) = Fh + X;
      end
      
      
      
     case 'whitened_posterior3'
      % Transform F to Z
      a = [covstruct.a(:); covstruct.s];
      V = linsolve_kron(@(x) feval(solvers{1,1}.linsolve, ...
                                   covstruct.LD_Cov{1}, ...
                                   x), ...
                        @(x) feval(solvers{2,1}.linsolve, ...
                                   covstruct.LD_Cov{2}, ...
                                   x), ...
                        Y) * sum(a)^2 / sum(a.^2);
      Z = zeros(size(F));
      for j=1:D_sum
        % Approximate posterior mean
        Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                         covstruct.K{2,j}, ...
                                         V);
% $$$         % Approximate whitening
% $$$         Z(:,:,j) = linsolve_kron(covstruct.L{1,j}, ...
% $$$                                  covstruct.L{2,j}, ...
% $$$                                  F(:,:,j) - Fh) / covstruct.a(j);
        % Approximate posterior whitening
        Z(:,:,j) = ...
            kronprod(covstruct.L_invpost{1,j}', ...
                     covstruct.L_invpost{2,j}', ...
                     linsolve_kron(@(x1) covstruct.L{1,j}\x1, ...
                                   @(x2) covstruct.L{2,j}\x2, ...
                                   F(:,:,j) - Fh)) / covstruct.a(j);
      end
      
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, Z, Imv);
      
      if ~isequal(covstruct, covstruct_new)
        covstruct = covstruct_new;
        rand_y = get_rand_y(covstruct);
      end
        
      % Transform Z back to F
      a = [covstruct.a(:); covstruct.s];
      V = linsolve_kron(@(x) feval(solvers{1,1}.linsolve, ...
                                   covstruct.LD_Cov{1}, ...
                                   x), ...
                        @(x) feval(solvers{2,1}.linsolve, ...
                                   covstruct.LD_Cov{2}, ...
                                   x), ...
                        Y) * sum(a)^2 / sum(a.^2);
      for j=1:D_sum
        % Approximate posterior mean (bias term)
        Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                         covstruct.K{2,j}, ...
                                         V);
        % Unwhiten Z
        X = linsolve_kron(covstruct.L_invpost{1,j}', ...
                          covstruct.L_invpost{2,j}', ...
                          Z(:,:,j));
        % Sum the unwhitened Z and the mean
        F(:,:,j) = Fh + kronprod(covstruct.L{1,j}, ...
                                 covstruct.L{2,j}, ...
                                 X) * covstruct.a(j);
      end
      
    end
    
    

    if debug
      disp('randfunc: Rand Y and F')
    end

% $$$     fprintf('Draw Y and F..\n');
% $$$     t = cputime();
% $$$     % Sample Y
% $$$     [Y, F] = rand_y(Y,Imv);
% $$$     %fprintf(' %f seconds\n', cputime()-t);
    
    % Process samples??
    if ~isempty(options.samplefunc)
      options.samplefunc(Y, F, covstruct.theta);
    end

    end

      function y = multiply_y(x)
      % Multiply matrix X by the full observation covariance matrix
      y = 0;
      for j=1:D_sum
        y = y + covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                            covstruct.K{2,j}, ...
                                            x);
      end
      y = y + s2W2.*x;
      end
  end
  
  function [covstruct, dcovstruct] = covstructfunc(theta, varargin)
  
% $$$   ind1 = 1+(1:n_theta1);
% $$$   ind2 = (1+n_theta1)+(1:n_theta2);
  
  thetas = cell(D_kron, D_sum);
  for j=1:D_sum
    for i=1:D_kron
      thetas{i,j} = theta(ind_theta{i,j});
    end
  end
  
  if debug
    disp('covstructfunc: start')
  end

  if nargin == 1
    
    covstruct.theta = theta;
    covstruct.a = theta(ind_a);
    covstruct.s = theta(ind_s);
    
    if nargout <= 1
      K = cell(D_kron, D_sum);
      for j=1:D_sum
        for i=1:D_kron
          K{i,j} = feval(covfuncs{i,j}, thetas{i,j});
        end
      end
    else
      K = cell(D_kron, D_sum);
      dK = cell(D_kron, D_sum);
      for j=1:D_sum
        for i=1:D_kron
          [K{i,j},dK{i,j}] = feval(covfuncs{i,j}, thetas{i,j});
        end
      end
    end
    
    covstruct.K = K;
    
    if nargout == 1
      L = cell(D_kron, D_sum);
      for j=1:D_sum
        for i=1:D_kron
          [L{i,j},p] = lchol(K{i,j});
          if p~=0
            L = NaN;
            break
          end
        end
        if ~iscell(L) && isnan(L)
          break;
        end
      end
      covstruct.L = L;
    else
      L = cell(D_kron, D_sum);
      dL = cell(D_kron, D_sum);
      for j=1:D_sum
        for i=1:D_kron
          [L{i,j},dL{i,j}] = lchol_grad(K{i,j},dK{i,j});
          if isnan(L{i,j})
            L = NaN;
            dL = NaN;
            break
          end
          %feval(solvers{i,j}.squareroot, covstruct.LD{i,j});
        end
        if ~iscell(L) && isnan(L)
          break
        end
      end
      covstruct.L = L;
      dcovstruct = covstruct;
      dcovstruct.dK = dK;
      dcovstruct.dL = dL;
    end      
    
  else
    
    if nargout > 1
      error('Not implemented for getting the gradient');
    end
    
    ind = varargin{1};
    covstruct_old = varargin{2};
    
    covstruct = covstruct_old;
    
    covstruct.theta = theta;
    covstruct.a = theta(ind_a);
    covstruct.s = theta(ind_s);
    
    for j=1:D_sum
      for i=1:D_kron
        if any(ismember(ind, ind_theta{i,j}))
          covstruct.K{i,j} = feval(covfuncs{i,j}, thetas{i,j});
          [covstruct.L{i,j},p] = lchol(covstruct.K{i,j});
          if p~=0
            covstruct.L{i,j} = NaN;
            break
          end
          %covstruct.LD{i,j} = feval(solvers{i,j}.decompose, covstruct.K{i,j});
        end
      end
    end
    
  end
  
% $$$   % Compute approximate covariance matrix
% $$$   covstruct.LD_Cov = cell(D_kron,1);
% $$$   for i=1:D_kron
% $$$     Cov = spdiag(covstruct.s^(2/D_kron) * w{i});
% $$$     for j=1:D_sum
% $$$       Cov = Cov + covstruct.a(j)^(2/D_kron) * covstruct.K{i,j};
% $$$     end
% $$$     covstruct.LD_Cov{i} = feval(solvers{i,1}.decompose, Cov);
% $$$   end
% $$$   
% $$$   % Compute approximate posterior covariances
% $$$   covstruct.L_invpost = cell(D_kron,D_sum);
% $$$   for i=1:D_kron
% $$$     invSigma = spdiag(1 ./ (covstruct.s^(2/D_kron)*w{i}));
% $$$     for j=1:D_sum
% $$$       Cov = covstruct.a(j)^(2/D_kron) ...
% $$$             * covstruct.L{i,j}' * invSigma * covstruct.L{i,j} ... 
% $$$             + speye(size(covstruct.L{i,j}));
% $$$       covstruct.L_invpost{i,j} = lchol(Cov);
% $$$     end
% $$$   end
  
  
  if debug
    disp('covstructfunc: end')
  end

  end
  
  function func = get_rand_y_pcg(covstruct)
  
  s2W2 = covstruct.s^2*W.^2;
  func = @rand_y_pcg;
    function [Y,F] = rand_y_pcg(Y, Imv)
    % Draw F from the prior
    F = zeros([size(Y),D_sum]);
    %size(F)
    for j=1:D_sum
      %size(covstruct.L{1,j})
      %size(covstruct.L{2,j})
      if isscalar(covstruct.L{1,j})
        foo = covstruct.L{1,j}
      end
      if isscalar(covstruct.L{2,j})
        bar = covstruct.L{2,j}
      end
      F(:,:,j) = covstruct.a(j) * gaussian_rand_kron(covstruct.L{1,j}, ...
                                                     covstruct.L{2,j});
    end
    % Draw Y (and subtract the observed matrix)
    Y0 = sum(F,3) + covstruct.s*W.*randn(size(Y)) - Y;
    % Solve CG once
    Z = conjgradmv(@multiply_y, ...
                   Y0, ...
                   ~Imv, ...
                   'maxiter', 1e3, ...
                   'tol', 1e-6, ...
                   'verbose', true, ...
                   'debug', false);
    % Transform prior samples to posterior samples
    for j=1:D_sum
      F(:,:,j) = F(:,:,j) - covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                                        covstruct.K{2,j}, ...
                                                        Z);
    end
    % Sample missing values
    sumF = sum(F,3);
    Y(Imv) = sumF(Imv) + ...
             covstruct.s*W(Imv) .* randn(sum(Imv(:)),1);
    
      function y = multiply_y(x)
      % Multiply matrix X by the full observation covariance matrix
      y = 0;
      for j=1:D_sum
        y = y + covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                            covstruct.K{2,j}, ...
                                            x);
      end
      y = y + s2W2.*x;
      end
    end
  end

  
  function func = get_loglikelihood_whitened_prior(covstruct)
  func = @loglikelihood_whitened_prior;
    function x = loglikelihood_whitened_prior(Y, Z, Imv)
    if ~iscell(covstruct.L) && isnan(covstruct.L)
      x = -inf;
      return
    end
    N_tot = sum(~Imv(:));
    ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(~Imv)));
    % Transform whitened variable V to the model variable F
    for j=1:D_sum
      if isnan(covstruct.L{1,j}(1))
        [j 1]
        x = -inf
        return
      end
      if isnan(covstruct.L{2,j}(1))
        [j 2]
        x = -inf
        return
      end
      Y = Y - covstruct.a(j) * kronprod(covstruct.L{1,j}, ...
                                        covstruct.L{2,j}, ...
                                        Z(:,:,j));
    end
    % Compute logpdf (ignore missing values)
    Y = Y ./ W;
    Y(Imv) = 0;
    x = gaussian_logpdf((Y(:)'*Y(:)) / (covstruct.s^2), ...
                        0, ...
                        0, ...
                        ldet_y, ...
                        N_tot);
    
    if isinf(x)
      x
      error('WHAAAT?!?! Logpdf is inf??')
    end
    end
  end
  
  function func = get_dloglikelihood_whitened_prior(dcovstruct)
  func = @dloglikelihood_whitened_prior;
    function dx = dloglikelihood_whitened_prior(Y, Z, Imv)
    if ~iscell(dcovstruct.L) && isnan(dcovstruct.L)
      dx = NaN;
      return
    end
    N_tot = sum(~Imv(:));
    ldet_y = 2*N_tot*log(dcovstruct.s) + 2*sum(log(W(~Imv)));
    % Transform whitened variable V to the model variable F
    sumF = zeros(size(Y));
    for j=1:D_sum
      sumF = sumF + dcovstruct.a(j) * kronprod(dcovstruct.L{1,j}, ...
                                         dcovstruct.L{2,j}, ...
                                         Z(:,:,j));
    end
    Y(Imv) = 0;
    sumF(Imv) = 0;
    Y = Y./W;
    sumF = sumF./W;
    % Initialize gradient
    dx = zeros(size(dcovstruct.theta));
    for j=1:D_sum
      % Gradient of scale
      V = kronprod(dcovstruct.L{1,j}, ...
                   dcovstruct.L{2,j}, ...
                   Z(:,:,j)) ./ W;
      dx(ind_a(j)) = ...
          gaussian_dlogpdf(0, ...
                           (Y(:)'*V(:)) / (dcovstruct.s^2), ...
                           2*(sumF(:)'*V(:)) / (dcovstruct.s^2), ...
                           0);
      % Gradient of covariance function parameters
      for i=1:length(ind_theta{1,j})
        V = dcovstruct.a(j) * kronprod(dcovstruct.dL{1,j}{i}, ...
                                       dcovstruct.L{2,j}, ...
                                       Z(:,:,j)) ./ W;
        dx(ind_theta{1,j}(i)) = ...
            gaussian_dlogpdf(0, ...
                             (Y(:)'*V(:)) / (dcovstruct.s^2), ...
                             2*(sumF(:)'*V(:)) / (dcovstruct.s^2), ...
                             0);
      end
      for i=1:length(ind_theta{2,j})
        V = dcovstruct.a(j) * kronprod(dcovstruct.L{1,j}, ...
                                       dcovstruct.dL{2,j}{i}, ...
                                       Z(:,:,j)) ./ W;
        dx(ind_theta{2,j}(i)) = ...
            gaussian_dlogpdf(0, ...
                             (Y(:)'*V(:)) / (dcovstruct.s^2), ...
                             2*(sumF(:)'*V(:)) / (dcovstruct.s^2), ...
                             0);
      end
    end
    % Gradient of noise term
    Y = Y-sumF;
    dx(ind_s) = gaussian_dlogpdf((Y(:)'*Y(:)) * (-2*dcovstruct.s^(-3)), ...
                                 0, ...
                                 0, ...
                                 2*N_tot/dcovstruct.s);
    end
  end
  
  function func = get_loglikelihood_whitened_posterior(covstruct)
  func = @loglikelihood_whitened_posterior;
    function x = loglikelihood_whitened_posterior(Y, Z, Imv)
    % Transform Z back to F
    a = [covstruct.a(:); covstruct.s];
    V = linsolve_kron(@(x) feval(solvers{1,1}.linsolve, ...
                                 covstruct.LD_Cov{1}, ...
                                 x), ...
                      @(x) feval(solvers{2,1}.linsolve, ...
                                 covstruct.LD_Cov{2}, ...
                                 x), ...
                      Y) * sum(a)^2 / sum(a.^2);
    x = gaussian_logpdf(Z(:)'*Z(:), ...
                        0, ...
                        0, ...
                        0, ...
                        numel(Z));
    for j=1:D_sum
      % Approximate posterior mean (bias term)
      Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                       covstruct.K{2,j}, ...
                                       V);
      % Unwhiten Z
      X = kronprod(covstruct.L{1,j}, ...
                   covstruct.L{2,j}, ...
                   Z(:,:,j)) * covstruct.a(j);
      % Sum the unwhitened Z and the bias
      F(:,:,j) = Fh + X;
      % Compute log p(Z)
      x = x + gaussian_logpdf(2*X(:)'*V(:) + Fh(:)'*V(:), ... 
                              0, ...
                              0, ...
                              0, ...
                              0);
      
    end
    
    % Compute log p(Y|Z), ignore missing values
    Y0 = (Y - sum(F,3)) ./ W;
% $$$     N_tot = sum(~Imv(:));
% $$$     ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(~Imv)));
% $$$     Y0(Imv) = 0;
    N_tot = numel(Y);
    ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(:)));
    x = x + gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.s^2), ...
                            0, ...
                            0, ...
                            ldet_y, ...
                            N_tot);
    
    end
  end
  
  function func = get_loglikelihood_whitened_posterior2(covstruct)
  s2W2 = covstruct.s^2*W.^2;
  func = @loglikelihood_whitened_posterior;
    function x = loglikelihood_whitened_posterior(Y, Z, Imv)
    % 
    V = conjgradmv(@multiply_y, ...
                   Y, ...
                   ~Imv, ...
                   'maxiter', 1e3, ...
                   'tol', 1e-6, ...
                   'verbose', true);
    
    x = gaussian_logpdf(Z(:)'*Z(:), ...
                        0, ...
                        0, ...
                        0, ...
                        numel(Z));
    for j=1:D_sum
      % Exact posterior mean (bias term)
      Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                       covstruct.K{2,j}, ...
                                       V);
      % Unwhiten Z
      X = kronprod(covstruct.L{1,j}, ...
                   covstruct.L{2,j}, ...
                   Z(:,:,j)) * covstruct.a(j);
      % Sum the unwhitened Z and the bias
      F(:,:,j) = Fh + X;
      % Compute log p(Z)
      x = x + gaussian_logpdf(2*X(:)'*V(:) + Fh(:)'*V(:), ... 
                              0, ...
                              0, ...
                              0, ...
                              0);
      
    end
    
    % Compute log p(Y|Z), ignore missing values
    Y0 = (Y - sum(F,3)) ./ W;
    N_tot = sum(~Imv(:));
    ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(~Imv)));
    Y0(Imv) = 0;
% $$$     N_tot = numel(Y);
% $$$     ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(:)));
    x = x + gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.s^2), ...
                            0, ...
                            0, ...
                            ldet_y, ...
                            N_tot);
    
      function y = multiply_y(x)
      % Multiply matrix X by the full observation covariance matrix
      y = 0;
      for j=1:D_sum
        y = y + covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                            covstruct.K{2,j}, ...
                                            x);
      end
      y = y + s2W2.*x;
      end
    end
  end
  
  
  function func = get_loglikelihood_whitened_posterior3(covstruct)
  func = @loglikelihood_whitened_posterior;
    function x = loglikelihood_whitened_posterior(Y, Z, Imv)
    
    % Approximate the posterior covariance with separable matrices
    
    % Variable for the approximate posterior mean
    a = [covstruct.a(:); covstruct.s];
    V = linsolve_kron(@(x) feval(solvers{1,1}.linsolve, ...
                                 covstruct.LD_Cov{1}, ...
                                 x), ...
                      @(x) feval(solvers{2,1}.linsolve, ...
                                 covstruct.LD_Cov{2}, ...
                                 x), ...
                      Y) * sum(a)^2 / sum(a.^2);
% $$$     x = gaussian_logpdf(Z(:)'*Z(:), ...
% $$$                         0, ...
% $$$                         0, ...
% $$$                         0, ...
% $$$                         numel(Z));
    x = 0;
    for j=1:D_sum
      % Approximate posterior mean (bias term)
      Fh = covstruct.a(j)^2 * kronprod(covstruct.K{1,j}, ...
                                       covstruct.K{2,j}, ...
                                       V);
      % Unwhiten Z
      X = linsolve_kron(covstruct.L_invpost{1,j}', ...
                        covstruct.L_invpost{2,j}', ...
                        Z(:,:,j));
      Xh = kronprod(covstruct.L{1,j}, ...
                    covstruct.L{2,j}, ...
                    X) * covstruct.a(j);
      % Compute log p(Z)
      ldet = N(2)*logdet_chol(covstruct.L_invpost{1,j}) ...
             + N(1)*logdet_chol(covstruct.L_invpost{2,j});
      x = x + gaussian_logpdf(Fh(:)'*V(:) + 2*V(:)'*Xh(:) + X(:)'*X(:), ... 
                              0, ...
                              0, ...
                              ldet, ...
                              numel(X));
      % Sum the unwhitened Z and the mean
      F(:,:,j) = Xh + Fh;
      
    end
    
    % Compute log p(Y|Z), ignore missing values
    Y0 = (Y - sum(F,3)) ./ W;
% $$$     N_tot = sum(~Imv(:));
% $$$     ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(~Imv)));
% $$$     Y0(Imv) = 0;
    N_tot = numel(Y);
    ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(:)));
    x = x + gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.s^2), ...
                            0, ...
                            0, ...
                            ldet_y, ...
                            N_tot);
    
    end
  end
  
end

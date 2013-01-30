
% K = SUM_d a_d^2*(K{1,d} * K{2,d}) + s^2*DIAG(W.^2)
%
% THETA = [A_1, THETA_11, THETA_21, A_2, THETA_12, THETA_22, ..., S]

function [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron6(covfuncs, solvers, logprior, dlogprior, varargin)


options = struct( ...
    'samplefunc',  [], ...
    'noise_scale', []);

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
if isempty(options.noise_scale)
  W = ones(N(2),N(1));
else
  W = options.noise_scale;
end

debug = false;

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

  covstruct = covstruct_init;
  
  func = @randfunc;
    function Y = randfunc(Y,Imv)
    
    % Draw a sample of the auxiliary variable
    Z = randn([N(end:-1:1), D_sum]);
  
    % Sample hyperparameters
    [~, covstruct] = rand_theta(Y, Z, Imv);

    % Process samples
    if ~isempty(options.samplefunc)
      feval(options.samplefunc, Y, F, covstruct.theta);
    end

    end

  end
  
  function [covstruct, dcovstruct] = covstructfunc(theta, varargin)
  
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
    
    covstruct.LD = cell(D_kron, D_sum);
    for j=1:D_sum
      for i=1:D_kron
        covstruct.LD{i,j} = feval(solvers{i,j}.decompose, covstruct.K{i,j});
      end
    end

    if nargout >= 2
      dcovstruct = covstruct;
      dcovstruct.dK = dK;
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
          covstruct.LD{i,j} = feval(solvers{i,j}.decompose, covstruct.K{i,j});
        end
      end
    end
    
  end
  
  % THIS STUFF IS FOR THE LIKELIHOOD EVALUATION AND SAMPLE DRAWING
  covstruct.L = cell(D_kron, D_sum);
  for j=1:D_sum
    for i=1:D_kron
      covstruct.L{i,j} = feval(solvers{i,j}.squareroot, covstruct.LD{i,j});
    end
  end
  
  if debug
    disp('covstructfunc: end')
  end

  end
  
  
  function func = get_loglikelihood(covstruct)
  L = cell(D_kron, D_sum);
  for j=1:D_sum
    for i=1:D_kron
      L{i,j} = feval(solvers{i,j}.squareroot, covstruct.LD{i,j});
    end
  end
  s2W2 = covstruct.s^2*W.^2;
  func = @loglikelihood_whitened_prior;
    function x = loglikelihood(Y, Z, Imv)
    
    % Initialize log-pdf value
    x = 0;
    
    % Draw a sample from the prior p(F|theta) using the auxiliary variable
    F_prior = zeros(size(F));
    for j=1:D_sum
      F_prior(:,:,j) = covstruct.a(j) * kronprod(covstruct.L{1,j}, ...
                                                 covstruct.L{2,j}, ...
                                                 Z(:,:,j));
    end
    % Draw a sample from the prior p(Y|F,theta) and subtract the observed Y
    Y0 = sum(F_prior,3) + covstruct.s * W .* randn(size(Y)) - Y;
    
    % Use conjugate gradient
    V = conjgradmv(@multiply_y, ...
                   Y0, ...
                   ~Imv, ...
                   'maxiter', 1e3, ...
                   'tol', 1e-6, ...
                   'verbose', true);
    
    % Compute cost from each p(F_j|theta) and transform prior samples to
    % posterior samples
    for j=1:D_sum
      % Compute Z - L'*V
      Z(:,:,j) = Z(:,:,j) - ...
          covstruct.a(j) * kronprod(covstruct.L{1,j}', ...
                                    covstruct.L{2,j}', ...
                                    V);
      % Compute log-determinant
      ldet = N(2) * feval(solvers{1,j}.logdet, covstruct.LD{1,j}) + ...
             N(1) * feval(solvers{2,j}.logdet, covstruct.LD{2,j});
      % Compute log-pdf (TODO: Maybe this doesn't need to be computed?)
      x = x + gaussian_logpdf(traceprod(Z,Z), ...
                              0 ...
                              0, ...
                              ldet, ...
                              prod(N));
      % Compute the posterior sample
      F_posterior(:,:,j) = covstruct.a(j) * kronprod(covstruct.L{1,j}, ...
                                                     covstruct.L{2,j}, ...
                                                     Z);
                                                     
    end
    
    % Compute cost from p(Y|F_j,theta)
    Y0 = (Y - sum(F_posterior,3)) ./ W;
    Y0(Imv) = 0;
    
    N_tot = sum(~Imv(:));
    ldet_y = 2*N_tot*log(covstruct.s) + 2*sum(log(W(~Imv)));
% $$$     % Transform whitened variable V to the model variable F
% $$$     F = zeros(size(Y));
% $$$     for j=1:D_sum
% $$$       F = F + covstruct.a(j) * kronprod(L{1,j},L{2,j},V(:,:,j));
% $$$     end
% $$$     % Compute logpdf (ignore missing values)
% $$$     Y0 = (Y - F) ./ W;
% $$$     Y0(Imv) = 0;
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
  
  
end

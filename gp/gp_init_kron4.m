
% k = a^2 * k1 * k2 + s^2
%
% K = a^2*kron(K1,K2) + s^2*eye(N1*N2)*diag(W(:).^2)
%
% THIS VERSION: The noise variances may be given fixed scales

function [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron4(covfunc1, solver1, covfunc2, solver2, logprior, dlogprior, ...
                  varargin)


options = struct( ...
    'samplefunc',  [], ...
    'likelihood',  'whitened_prior', ...
    'noise_scale', [], ...
    'rand_y',      'pcg'); % gibbs / pcg

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

[n_theta1, N1, M1] = covfunc1();
[n_theta2, N2, M2] = covfunc2();

if isempty(options.noise_scale)
  W = ones(N2,N1);
else
  W = reshape(options.noise_scale, N2, N1);
end

debug = false;

% Posterior density function
switch options.likelihood
 case 'conditioned'
  get_loglikelihood = @get_loglikelihood_conditioned;
  get_dloglikelihood = @get_dloglikelihood_conditioned;
 case 'whitened_prior'
  get_loglikelihood = @get_loglikelihood_whitened_prior;
  get_dloglikelihood = @get_dloglikelihood_whitened_prior;
% $$$  case 'whitened_approximate'
% $$$   get_loglikelihood = @get_loglikelihood_whitened_approximate;
% $$$   get_dloglikelihood = @get_dloglikelihood_whitened_approximate;
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
  
    if isempty(F)
      % Use some better initialization here..
      fprintf('Draw Y and F..\n');
      t = cputime();
      [Y, F] = rand_y(Y,Imv);
      %fprintf(' %f seconds\n', cputime()-t);
    end
    
    fprintf('Draw theta..\n');
    t = cputime();

    switch options.likelihood
      
     case 'conditioned'
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, F);
      
     case 'whitened_prior'
      % Prior whitening, transform to V..
      U1 = solver1.squareroot(covstruct.L1)';
      U2 = solver2.squareroot(covstruct.L2)';
      V = kronprod(U1, ...
                   U2, ...
                   linsolve_kron(@(x1) solver1.linsolve(covstruct.L1,x1), ...
                                 @(x2) solver2.linsolve(covstruct.L2,x2), ...
                                 F)) / covstruct.a;
      
      disp('Whitened prior')
      
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, V);
      
% $$$       % .. and transform back to F
% $$$       F = covstruct.a * kronprod(U1', U2', V);
% $$$       Y(Imv) = F(Imv) + covstruct.s*W(Imv).*randn(size(F(Imv)));

    end
    
    if ~isequal(covstruct, covstruct_new)
      covstruct = covstruct_new;
      rand_y = get_rand_y(covstruct);
      % If whitened likelihood, also update F..
    end
    

    if debug
      disp('randfunc: Rand Y and F')
    end

    fprintf('Draw Y and F..\n');
    t = cputime();
    % Sample Y
    [Y, F] = rand_y(Y,Imv);
    %fprintf(' %f seconds\n', cputime()-t);
    
    % Process samples??
    if ~isempty(options.samplefunc)
      options.samplefunc(Y, F, covstruct.theta);
    end

    end

  end
  
  function [covstruct, dcovstruct] = covstructfunc(theta, varargin)
  
  ind1 = 1+(1:n_theta1);
  ind2 = (1+n_theta1)+(1:n_theta2);
  
  theta1 = theta(ind1);
  theta2 = theta(ind2);
  
  if debug
    disp('covstructfunc: start')
  end

  if nargin == 1
    
    covstruct.theta = theta;
    covstruct.a = theta(1);
    covstruct.s = theta(end);
    
    if nargout <= 1
      K1 = covfunc1(theta1);
      K2 = covfunc2(theta2);
    else
      [K1, dK1] = covfunc1(theta1);
      [K2, dK2] = covfunc2(theta2);
    end
    
    covstruct.K1 = K1;
    covstruct.K2 = K2;

    covstruct.L1 = solver1.decompose(covstruct.K1);
    covstruct.L2 = solver2.decompose(covstruct.K2);
    
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
    
    covstruct.theta = theta;
    covstruct.a = theta(1);
    covstruct.s = theta(end);
    
    if any(ismember(ind, ind1))
      covstruct.K1 = covfunc1(theta1);
      covstruct.L1 = solver1.decompose(covstruct.K1);
    end
    if any(ismember(ind, ind2))
      covstruct.K2 = covfunc2(theta2);
      covstruct.L2 = solver2.decompose(covstruct.K2);
    end
    
  end

  if debug
    disp('covstructfunc: end')
  end

  end
  
  function func = get_rand_y_pcg(covstruct)
  L1 = solver1.squareroot(covstruct.L1);
  L2 = solver2.squareroot(covstruct.L2);
  multiply_f = @(x) covstruct.a^2*kronprod(covstruct.K1, covstruct.K2, x);
  multiply_y = @(x) (multiply_f(x) + covstruct.s^2*W.^2.*x);
  func = @rand_y_pcg;
    function [Y,F] = rand_y_pcg(Y, Imv)
    F = covstruct.a * gaussian_rand_kron(L1,L2);
    F = gaussian_rand_conjgradmv(multiply_f, ...
                                 multiply_y, ...
                                 F, ...
                                 F + covstruct.s*W.*randn(size(Y)) - Y, ...
                                 ~Imv, ...
                                 'maxiter', 1e3, ...
                                 'tol', 1e-6, ...
                                 'verbose', true);
    Y(Imv) = F(Imv) + covstruct.s*W(Imv).*randn(size(Y(Imv)));
% $$$     if ~isreal(F)
% $$$       covstruct.L1.LD
% $$$       full(L1)
% $$$       error('Hmm.. F not real..');
% $$$     end
    end
% $$$   L1 = solver1.squareroot(covstruct.L1);
% $$$   L2 = solver2.squareroot(covstruct.L2);
% $$$   multiply_f = @(x) covstruct.a^2*reshape(kronprod(covstruct.K1, ...
% $$$                                                    covstruct.K2, ...
% $$$                                                    reshape(x,[N2 N1])), ...
% $$$                                           [N1*N2,1]);
% $$$   w2 = W(:).^2;
% $$$   multiply = @(x) (multiply_f(x) + covstruct.s^2*w2.*x);
% $$$   M = diag(covstruct.K2)*diag(covstruct.K1)' + W.^2;
% $$$   preconditioner = @(x) x./M(:);
% $$$   func = @rand_y_pcg;
% $$$     function [Y,F] = rand_y_pcg(Y, Imv)
% $$$     S = speye(numel(Y));
% $$$     S(Imv(:),:) = [];
% $$$     Z_f = covstruct.a * gaussian_rand_kron(L1,L2);
% $$$     Z_noise = covstruct.s * W(:) .* randn(numel(Y),1);
% $$$     [F] = gaussian_rand_pcg(Y(:), ...
% $$$                                     multiply, ...
% $$$                                     multiply_f, ...
% $$$                                     Z_f(:), ...
% $$$                                     Z_noise(:), ...
% $$$                                     preconditioner, ...
% $$$                                     S, ...
% $$$                                     'maxiter', 1e2, ...
% $$$                                     'tol', 1e-3);
% $$$     F = reshape(F, size(Y));
% $$$     %F_mean = reshape(F_mean, size(Y));
% $$$     Y(Imv) = F(Imv) + covstruct.s * W(Imv) .* randn(size(Y(Imv)));
% $$$     
% $$$     end
  end

  
  function func = get_loglikelihood_whitened_prior(covstruct)
  L1 = solver1.squareroot(covstruct.L1);
  L2 = solver2.squareroot(covstruct.L2);
  ldet_f = 0; %N2*ldet1 + N1*ldet2 + 2*N2*N1*log(covstruct.a);
  N = N1*N2;
  ldet_y = 2*N*log(covstruct.s) + 2*sum(log(W(:)));
  func = @loglikelihood_whitened_prior;
    function x = loglikelihood_whitened_prior(Y, V)
    Y0 = (Y - kronprod(L1,L2,V)*covstruct.a) ./ W;
    % Note that the first term doesn't depend on covstruct, i.e., the
    % hyperparameters..
    % Why don't we just ignore missing values?
    x = gaussian_logpdf((V(:)'*V(:)), ...
                        0, ...
                        0, ...
                        ldet_f, ...
                        N) + ...
        gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.s^2), ...
                        0, ...
                        0, ...
                        ldet_y, ...
                        N);
    end
  end
  
  
end

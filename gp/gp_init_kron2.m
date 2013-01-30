
% k = a^2 * k1 * k2 + s^2
%
% K = a^2*kron(K1,K2) + s^2*eye(N1*N2)

function [get_logpdf, get_dlogpdf, get_rand_model, func_theta] = ...
    gp_init_kron2(covfunc1, solver1, covfunc2, solver2, logprior, dlogprior, ...
                  varargin)


options = struct( ...
    'samplefunc', [], ...
    'likelihood', 'whitened_approximate', ...
    'rand_y',  'gibbs'); % gibbs / pcg

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

[n_theta1, N1, M1] = covfunc1();
[n_theta2, N2, M2] = covfunc2();

debug = false;

% Allocate memory
%TMP = zeros(N2, N1);

% Posterior density function
switch options.likelihood
 case 'conditioned'
  get_loglikelihood = @get_loglikelihood_conditioned;
  get_dloglikelihood = @get_dloglikelihood_conditioned;
 case 'whitened_prior'
  get_loglikelihood = @get_loglikelihood_whitened_prior;
  get_dloglikelihood = @get_dloglikelihood_whitened_prior;
 case 'whitened_approximate'
  get_loglikelihood = @get_loglikelihood_whitened_approximate;
  get_dloglikelihood = @get_dloglikelihood_whitened_approximate;
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
      % Prior whitening
      U1 = solver1.squareroot(covstruct.L1)';
      U2 = solver2.squareroot(covstruct.L2)';
      V = kronprod(U1, ...
                   U2, ...
                   linsolve_kron(@(x1) solver1.linsolve(covstruct.L1,x1), ...
                                 @(x2) solver2.linsolve(covstruct.L2,x2), ...
                                 F)) / covstruct.a;
% $$$       L1 = solver1.squareroot(covstruct.L1);
% $$$       L2 = solver2.squareroot(covstruct.L2);
% $$$       V = linsolve_kron(@(x) linsolve_tril(L1,x), ...
% $$$                         @(x) linsolve_tril(L2,x), ...
% $$$                         F) / covstruct.a;
      
      disp('Whitened prior')
      
      % Sample hyperparameters
      [~, covstruct_new] = rand_theta(Y, V);

     case 'whitened_approximate'
      
      if debug
        disp('randfunc: Evaluate V')
      end
      
      linsolve_L = ...
          @(LD1,LD2,X) kronprod(solver1.squareroot(LD1)', ...
                                solver2.squareroot(LD2)', ...
                                linsolve_kron(@(x1) solver1.linsolve(LD1,x1), ...
                                              @(x2) solver2.linsolve(LD2,x2), ...
                                              X));
      
      U_Cov1 = solver1.squareroot(covstruct.L_Cov1)';
      U_Cov2 = solver2.squareroot(covstruct.L_Cov2)';
      %L1 = solver1.squareroot(covstruct.L1);
      %L2 = solver2.squareroot(covstruct.L2);
      
      linsolve_cov1 = @(x) solver1.linsolve(covstruct.L_Cov1, x);
      linsolve_cov2 = @(x) solver2.linsolve(covstruct.L_Cov2, x);

      V = kronprod(U_Cov1, ...
                   U_Cov2, ...
                   linsolve_L(covstruct.L1, ...
                              covstruct.L2, ...
                              F - covstruct.a^2*kronprod(covstruct.K1, ...
                                                        covstruct.K2, ...
                                                        linsolve_kron(linsolve_cov1, ...
                                                        linsolve_cov2, ...
                                                        Y) ) ...
                              ) ...
                   ) / (covstruct.s*covstruct.a);
% $$$       V = kronprod(U_Cov1, ...
% $$$                    U_Cov2, ...
% $$$                    linsolve_kron(L1, ...
% $$$                                  L2, ...
% $$$                                  F - covstruct.a^2*kronprod( covstruct.K1, ...
% $$$                                                covstruct.K2, ...
% $$$                                                linsolve_kron(linsolve_cov1, ...
% $$$                                                         linsolve_cov2, ...
% $$$                                                         Y) ) ...
% $$$                                  ) ...
% $$$                    ) / (covstruct.s*covstruct.a);
      
      disp('Whitened approximate')
      
      % Sample hyperparameters
      if debug
        disp('randfunc: Rand theta')
      end
      [~, covstruct_new] = rand_theta(Y, V);

      if debug
        disp('randfunc: Evaluate F')
      end

      if ~isequal(covstruct, covstruct_new)
        covstruct = covstruct_new;
        rand_y = get_rand_y(covstruct);
        % If whitened likelihood, also update F..
      end
      
      linsolve_U = ...
          @(LD1,LD2,X) linsolve_kron(@(x1) solver1.linsolve(LD1,x1), ...
                                     @(x2) solver2.linsolve(LD2,x2), ...
                                     kronprod(solver1.squareroot(LD1),...
                                              solver2.squareroot(LD2), ...
                                              X));

      %U_Cov1 = solver1.squareroot(covstruct.L_Cov1)';
      %U_Cov2 = solver2.squareroot(covstruct.L_Cov2)';
      L1 = solver1.squareroot(covstruct.L1);
      L2 = solver2.squareroot(covstruct.L2);
      
      linsolve_cov1 = @(x) solver1.linsolve(covstruct.L_Cov1, x);
      linsolve_cov2 = @(x) solver2.linsolve(covstruct.L_Cov2, x);

      % Back to F
      F = covstruct.a^2*kronprod(covstruct.K1, ...
                                 covstruct.K2, ...
                                 linsolve_kron(linsolve_cov1, ...
                                               linsolve_cov2, ...
                                               Y)) + ...
          (covstruct.a*covstruct.s) * kronprod(L1, ...
                                               L2, ...
                                               linsolve_U(covstruct.L_Cov1, ...
                                                        covstruct.L_Cov2, ...
                                                        V));
      % ... and Y
      Y(Imv) = F(Imv) + covstruct.s*randn(size(F(Imv)));
% $$$       Y(Imv) = F(Imv) + covstruct.a*covstruct.s*randn(size(F(Imv)));
      
    end
    
    if ~isequal(covstruct, covstruct_new)
      covstruct = covstruct_new;
      rand_y = get_rand_y(covstruct);
      % If whitened likelihood, also update F..
    end
    

    %fprintf(' %f seconds\n', cputime()-t);
    
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
    
    %%%%
    Cov1 = covstruct.a*covstruct.K1 + covstruct.s*speye(N1);
    Cov2 = covstruct.a*covstruct.K2 + covstruct.s*speye(N2);
% $$$     Cov1 = covstruct.K1 + covstruct.s*speye(N1);
% $$$     Cov2 = covstruct.K2 + covstruct.s*speye(N2);
    covstruct.L_Cov1 = solver1.decompose(Cov1);
    covstruct.L_Cov2 = solver2.decompose(Cov2);
    %%%%
    
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
    
    %%%%
    if any(ind==1) || any(ismember(ind,ind1)) || any(ind==length(theta))
      Cov1 = covstruct.a*covstruct.K1 + covstruct.s*speye(N1);
      covstruct.L_Cov1 = solver1.decompose(Cov1);
    end
    if any(ind==1) || any(ismember(ind,ind2)) || any(ind==length(theta))
      Cov2 = covstruct.a*covstruct.K2 + covstruct.s*speye(N2);
      covstruct.L_Cov2 = solver2.decompose(Cov2);
    end
% $$$     if any(ismember(ind,ind1)) || any(ind==length(theta))
% $$$       Cov1 = covstruct.K1 + covstruct.s*speye(N1);
% $$$       covstruct.L_Cov1 = solver1.decompose(Cov1);
% $$$     end
% $$$     if any(ismember(ind,ind2)) || any(ind==length(theta))
% $$$       Cov2 = covstruct.K2 + covstruct.s*speye(N2);
% $$$       covstruct.L_Cov2 = solver2.decompose(Cov2);
% $$$     end
    %%%%
  
  end

  if debug
    disp('covstructfunc: end')
  end

  end
  
  function func = get_rand_y_pcg(covstruct)
  L1 = solver1.squareroot(covstruct.L1);
  L2 = solver2.squareroot(covstruct.L2);
  multiply_f = @(x) covstruct.a^2*kronprod(covstruct.K1, covstruct.K2, x);
% $$$   multiply_f = @(x) covstruct.a^2*reshape(kronprod(covstruct.K1, ...
% $$$                                                    covstruct.K2, ...
% $$$                                                    reshape(x,[N2 N1])), ...
% $$$                                           [N1*N2,1]);
  multiply_y = @(x) (multiply_f(x) + covstruct.s^2*x);
  func = @rand_y_pcg;
    function [Y,F] = rand_y_pcg(Y, Imv)
    %S = speye(numel(Y));
    %S(Imv(:),:) = [];
    F = covstruct.a * gaussian_rand_kron(L1,L2);
    %Z_noise = covstruct.s * randn(size(Y));
% $$$     disp('keyboard in rand_y_pcg');
% $$$     keyboard()
    F = gaussian_rand_conjgradmv(multiply_f, ...
                                 multiply_y, ...
                                 F, ...
                                 F + covstruct.s*randn(size(Y)) - Y, ...
                                 ~Imv, ...
                                 'maxiter', 1e3, ...
                                 'tol', 1e-6, ...
                                 'verbose', true);
% $$$     F = gaussian_rand_conjgradmv(Y, ...
% $$$                                  multiply, ...
% $$$                                  multiply_f, ...
% $$$                                  F, ...
% $$$                                  Z_y, ...
% $$$                                  ~Imv, ...
% $$$                                  'maxiter', 1e2, ...
% $$$                                  'tol', 1e-6, ...
% $$$                                  'verbose', true);
    %keyboard
    %F = reshape(F, size(Y));
    Y(Imv) = F(Imv) + covstruct.s*randn(size(Y(Imv)));
% $$$       function x = multiply_f_pcg(x)
% $$$       A = zeros(size(Y));
% $$$       A(~Imv) = x;
% $$$       A = multiply(A);
% $$$       x = A;
% $$$       end
% $$$       function x = multiply_pcg(x)
% $$$       A = zeros(size(Y));
% $$$       A(~Imv) = x;
% $$$       A = multiply(A);
% $$$       x = A(~Imv);
% $$$       end
    end
  end

  function func = get_rand_y_gibbs(covstruct)

  L1 = solver1.squareroot(covstruct.L1);
  L2 = solver2.squareroot(covstruct.L2);
  
  % Noiseless part of the covariance
  multiply_f = @(X) covstruct.a^2*kronprod(covstruct.K1,covstruct.K2,X);
  
  % Linear system for CG
  multiply = @(X) (multiply_f(X)+covstruct.s^2*X);
% $$$   multiply = @(X) (multiply_f(X)+covstruct.a^2*covstruct.s^2*X);
  multiply_pcg = @(X) reshape(multiply(reshape(X, [N2 N1])), [N2*N1,1]);

  % Preconditioner for CG
  linsolve1_precond = @(x) solver1.linsolve(covstruct.L_Cov1, x);
  linsolve2_precond = @(x) solver2.linsolve(covstruct.L_Cov2, x);
  linsolve_precond = @(X) linsolve_kron(linsolve1_precond, ...
                                        linsolve2_precond, ...
                                        X);
% $$$   linsolve_precond = @(X) linsolve_kron(linsolve1_precond, ...
% $$$                                         linsolve2_precond, ...
% $$$                                         X) / covstruct.a^2;
  linsolve_precond_pcg = ...
      @(X) reshape(linsolve_precond(reshape(X, [N2 N1])), ...
                   [N2*N1,1]);
  
  % CG solver
  linsolve = @(X) reshape(pcg(multiply_pcg, ...
                              X(:), ...
                              1e-3, ...
                              100,  ...
                              linsolve_precond_pcg, ...
                              []), ...
                          [N2, N1] );
  

  func = @rand_y_gibbs;

    function [Y, F] = rand_y_gibbs(Y, Imv)
    % Sample latent function values
    if debug
      disp('rand_y_gibbs: 1')
    end
    TMP = covstruct.a*gaussian_rand_kron(L1,L2);
    if debug
      disp('rand_y_gibbs: 2')
    end
    F = TMP - multiply_f( linsolve( TMP + covstruct.s*randn(N2,N1) - Y ) );
% $$$     F = TMP - multiply_f( linsolve( TMP + covstruct.a*covstruct.s*randn(N2,N1) - Y ) );
    % Sample missing values
    if debug
      disp('rand_y_gibbs: 3')
    end
    Y(Imv) = F(Imv) + covstruct.s*randn(size(Y(Imv)));
% $$$     Y(Imv) = F(Imv) + covstruct.a*covstruct.s*randn(size(Y(Imv)));
    % You could evaluate:
    % rand(F), mean(F), mean(Y), 
    end
  
  end
  
  function func = get_loglikelihood_whitened_approximate(covstruct)

  %U_Cov1 = solver1.squareroot(covstruct.L_Cov1)';
  %U_Cov2 = solver2.squareroot(covstruct.L_Cov2)';
  L1 = solver1.squareroot(covstruct.L1);
  L2 = solver2.squareroot(covstruct.L2);
  
  linsolve_cov1 = @(x) solver1.linsolve(covstruct.L_Cov1, x);
  linsolve_cov2 = @(x) solver2.linsolve(covstruct.L_Cov2, x);
  
  linsolve = @(X) linsolve_kron(@(x1) solver1.linsolve(covstruct.L1, x1), ...
                                @(x2) solver2.linsolve(covstruct.L2, x2), ...
                                X);
  ldet1 = solver1.logdet(covstruct.L1);
  ldet2 = solver2.logdet(covstruct.L2);
  ldet_f = N2*ldet1 + N1*ldet2 + 2*N2*N1*log(covstruct.a);
  N = N1*N2;
  ldet_y = 2*N*log(covstruct.s);
% $$$   ldet_y = 2*N*log(covstruct.s) + 2*N*log(covstruct.a);

  ldet_jacob = - 2*N*log(covstruct.s) - 2*N*log(covstruct.a) - N2*ldet1 - N1*ldet2...
      + N2*solver1.logdet(covstruct.L_Cov1) + N1*solver2.logdet(covstruct.L_Cov2);
  
  linsolve_U = ...
      @(LD1,LD2,X) linsolve_kron(@(x1) solver1.linsolve(LD1,x1), ...
                                 @(x2) solver2.linsolve(LD2,x2), ...
                                 kronprod(solver1.squareroot(LD1),...
                                          solver2.squareroot(LD2), ...
                                          X));

  func = @loglikelihood_whitened_approximate;
    function x = loglikelihood_whitened_approximate(Y, V)
    if debug
      disp('loglikelihood_whitened_approximate: 1')
    end
    F = covstruct.a^2*kronprod(covstruct.K1, ...
                               covstruct.K2, ...
                               linsolve_kron(linsolve_cov1, ...
                                             linsolve_cov2, ...
                                             Y)) + ...
        covstruct.a*covstruct.s * kronprod(L1, ...
                                           L2, ...
                                           linsolve_U(covstruct.L_Cov1, ...
                                                      covstruct.L_Cov2, ...
                                                      V));
    if debug
      disp('loglikelihood_whitened_approximate: 2')
    end
    TMP = Y-F;
    if debug
      disp('loglikelihood_whitened_approximate: 3')
    end
    x = gaussian_logpdf(traceprod(F, linsolve(F)) / covstruct.a^2, ...
                        0, ...
                        0, ...
                        ldet_f + ldet_jacob, ...
                        N) + ...
        gaussian_logpdf((TMP(:)'*TMP(:)) / (covstruct.s^2), ...
                        0, ...
                        0, ...
                        ldet_y, ...
                        N);
% $$$         gaussian_logpdf((TMP(:)'*TMP(:)) / (covstruct.a^2*covstruct.s^2), ...
% $$$                         0, ...
% $$$                         0, ...
% $$$                         ldet_y, ...
% $$$                         N);
    end
  end
  
  function func = get_loglikelihood_whitened_prior(covstruct)
  L1 = solver1.squareroot(covstruct.L1);
  L2 = solver2.squareroot(covstruct.L2);
  ldet_f = 0; %N2*ldet1 + N1*ldet2 + 2*N2*N1*log(covstruct.a);
  N = N1*N2;
  ldet_y = 2*N*log(covstruct.s);
% $$$   ldet_y = 2*N*log(covstruct.s) + 2*N*log(covstruct.a);
  func = @loglikelihood_whitened_prior;
    function x = loglikelihood_whitened_prior(Y, V)
    Y0 = (Y - kronprod(L1,L2,V)*covstruct.a);
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
% $$$         gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.a^2*covstruct.s^2), ...
% $$$                         0, ...
% $$$                         0, ...
% $$$                         ldet_y, ...
% $$$                         N);
    end
  end
  
  function func = get_loglikelihood_conditioned(covstruct)
  linsolve = @(X) linsolve_kron(@(x1) solver1.linsolve(covstruct.L1, x1), ...
                                @(x2) solver2.linsolve(covstruct.L2, x2), ...
                                X);
  ldet1 = solver1.logdet(covstruct.L1);
  ldet2 = solver2.logdet(covstruct.L2);
  ldet_f = N2*ldet1 + N1*ldet2 + 2*N2*N1*log(covstruct.a);
  N = N1*N2;
  ldet_y = 2*N*log(covstruct.s);
% $$$   ldet_y = 2*N*log(covstruct.s) + 2*N*log(covstruct.a);
  func = @loglikelihood_conditioned;
    function x = loglikelihood_conditioned(Y, F)
    Y0 = (Y-F);
    x = gaussian_logpdf(traceprod(F, linsolve(F)) / covstruct.a^2, ...
                        0, ...
                        0, ...
                        ldet_f, ...
                        N) + ...
        gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.s^2), ...
                        0, ...
                        0, ...
                        ldet_y, ...
                        N);
% $$$         gaussian_logpdf((Y0(:)'*Y0(:)) / (covstruct.a^2*covstruct.s^2), ...
% $$$                         0, ...
% $$$                         0, ...
% $$$                         ldet_y, ...
% $$$                         N);
    end
  end
  
  function func = get_dloglikelihood_conditioned(dcovstruct)
  warning('Check the gradients - scale covstruct.a is not correct here')
  linsolve = @(X) linsolve_kron(@(x1) solver1.linsolve(dcovstruct.L1, x1), ...
                                @(x2) solver2.linsolve(dcovstruct.L2, x2), ...
                                X);
  invK1 = solver1.inv(dcovstruct.L1);
  invK2 = solver2.inv(dcovstruct.L2);
  N = N1*N2;
  func = @dloglikelihood_conditioned;
    function dx = dloglikelihood_conditioned(Y, F)
    Z_f = linsolve(F);
    Y0 = (Y-F);
    n_theta = 1 + n_theta1 + n_theta2 + 1;
    dx = zeros(n_theta,1);
    % Gradient for scale parameter
    dx(1) = - 0.5 * traceprod(F, Z_f) * (-2)*dcovstruct.a^(-3) ...
            - N / dcovstruct.a;
    % Gradient for K1 hyperparameters
    for n=1:n_theta1
      n1 = 1 + n;
      dx(n1) = 0.5 * traceprod(Z_f,kronprod(dcovstruct.dK1{n}, ...
                                            dcovstruct.K2, ...
                                            Z_f)) / dcovstruct.a^2 - ...
               0.5 * N2 * traceprod(invK1, ...
                                    dcovstruct.dK1{n});
    end
    % Gradient for K2 hyperparameters
    for n=1:n_theta2
      n2 = 1 + n_theta1 + n;
      dx(n2) = 0.5 * traceprod(Z_f,kronprod(dcovstruct.K1, ...
                                            dcovstruct.dK2{n}, ...
                                            Z_f)) / dcovstruct.a^2 - ...
               0.5 * N1 * traceprod(invK2, ...
                                    dcovstruct.dK2{n});
    end
    % Gradient for noise parameter
    dx(end) = - 0.5 * traceprod(Y0,Y0) * (-2)*dcovstruct.s.^(-3) - ...
              N / dcovstruct.s;
    end
  end
  
end

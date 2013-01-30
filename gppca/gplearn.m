%
% [f, Covf, pseudoX, logtheta, Kpp, Kxp, Lp] = 
%     gplearn(y, V, X, pseudoX, kX, logtheta, ...)
%
% V is inverse covariance (i.e. precision) matrix for observation noise -
% sparse format is recommended.
%
% In practice, inv(V) should be diagonal. However, inv(V) can be more
% complex sparse matrix as long as it has zeros where kX is non-zero. That
% is, observation noise can be correlated if the corresponding latent
% function values are independent a priori. This more general property is
% used in PCA model, where inv(V) and kX are block-diagonal (although rows
% and columns may be in mixed order!) but the blocks "overlap" only on the
% diagonal, that is, they both have non-zero elements at the same time only
% on the diagonal.
%
% y = inv(V)*c, that is, c = V*y
% N( inv(V)*c | 0, inv(V) + Kxp inv(Kpp) Kpx )

% $$$ function [f, Covf, pseudoX, logtheta, Kpp, Kxp, Lp] = gplearn(y, V, ...
% $$$                                                   X, pseudoX, kX, logtheta, ...
% $$$                                                   varargin)
function [f, Covf, Xp, logtheta, Kpp, Kxp, Lp] = gplearn(logtheta, covfunc, ...
                                                  X, y, V, Xp, varargin)

% function [] = gp_learn(covfunc, theta, y, V)

warning('Deprecated. Use gp_learn or gp_learn_pseudo');


if ~iscell(covfunc)
  covfunc = {covfunc};
end

opts = struct( ...
    'maxsearch', 10, ...
    'vy', [], ...
    'cholv', [], ...
    'checkgrad', false, ...
    'updatepseudo', true, ...
    'priorlogtheta', []);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

m = cols(Xp); % number of pseudoinputs
d = rows(X); % dimensionality of the input space

if isempty(opts.cholv)
  LV = chol(V, 'lower');
%  [LV,V] = safechol(V, 1e-16, 'lower');
else
  LV = opts.cholv;
end

if isempty(opts.vy)
  c = V*y;
else
  c = opts.vy;
end

if opts.maxsearch ~= 0
  if opts.updatepseudo && ~isempty(Xp)
    % Learn pseudo-inputs and hyperparameters
    x0 = [Xp(:); logtheta(:)];
    func = @(x) bound2cost(x, c, V, LV, X, covfunc, m, d, true);
    fiter = inf;
    while ~any(~isinf(fiter))
      [x, fiter] = minimize(x0, func, opts.maxsearch);
      x0 = x0 - 1; % adhoc trying to help ill conditioned matrices..
    end
    Xp = reshape(x(1:(m*d)), [d,m]);
    logtheta = x((m*d+1):end);
  else
    % Learn only hyperparameters
    x0 = logtheta(:);
    func = @(x) bound2cost(x, c, V, LV, X, covfunc, [], d, false, Xp);
    fiter = inf;
    while ~any(~isinf(fiter))
      [logtheta, fiter] = minimize(x0, func, opts.maxsearch);
      x0 = x0 - 1; % adhoc trying to help ill conditioned matrices..
    end
  end
  if opts.checkgrad
    mycheckgrad(func, x0, 1e-8)
  end
end

% THE FOLLOWING WOULD COME AS A SIDE PRODUCT FROM THE OPTIMIZATION BUT HOW
% TO GET THEM??

if isempty(Xp)
  % Full case (no pseudo-inputs)
  
  n = cols(X);
  In = diag(sparse(ones(n,1)));

  K = feval(covfunc{:}, logtheta, X, X);
% $$$   norm_K = normest(K);
% $$$   K = K + 1e-14*norm_K*In; % SET PROPERLY
%  K = K + 1e-6*In;

  % DEBUG: REGULARIZE
  K = regularize(K);
  Lp = chol(K, 'lower');
%  [Lp,K] = safechol(K, 1e-10, 'lower');
  %K = Lp*Lp';
  Kpp = K;
  Kxp = [];
  
  % We assume that V is DIAGONAL
  V = sparse(V);
  L = chol(inv(V) + K, 'lower');
%  L = safechol(inv(V) + K, 1e-16, 'lower');
  
  % Posterior mean and covariance
  if false
    % TODO: NOTE, THAT DIAGONAL COVARIANCE APPROXIMATION IS A BIT
    % ADHOC!!!!! CAUSES LOGLIKELIHOOD TO BE INCORRECT!!!
    f = K * c - K * solve_triu(L', solve_tril(L, K*c));
    Covf = spalloc(n,n,n);
    for i=1:n
      r = solve_tril(L,K(:,i));
      Covf(i,i) = K(i,i) - r'*r;
    end
  else
    invLK = full(solve_tril(L,K));
    Covf = K - invLK' * invLK;
    f = Covf * c;
  end
  
else
  % Pseudo input case

  Im = diag(sparse(ones(m,1)));
  
  % GP-prior covariances
  Kpp = feval(covfunc{:}, logtheta, Xp, Xp);
% $$$   norm_Kpp = normest(Kpp);
% $$$   Kpp = Kpp + 1e-14*norm_Kpp*Im; % SET PROPERLY
%  Kpp = Kpp + 1e-3*Im; % SET PROPERLY
  Kxp = feval(covfunc{:}, logtheta, X, Xp);
  
  % DEBUG: REGULARIZE
  %Kpp = regularize(Kpp);

  % Helpful variables
  reg = 1;
  while 1
    try
      Kpp = regularize(Kpp, reg);
      Lp = chol(Kpp, 'lower');
      break
    catch
      warning('Ill conditioned covariance matrix in learning phase!!');
      reg = reg * 10;
    end
  end
  %  [Lp, Kpp] = safechol(Kpp, 1e-10, 'lower');
  %Kpp = Lp*Lp';

  %szLp = size(Lp)
  %szKxp = size(Kxp)
  %szLV = size(LV)
  R = solve_tril(Lp, Kxp' * LV);
  S = Im + R*R';
  LS = chol(S, 'lower');
%  LS = safechol(S, 1e-16, 'lower');
  z = LV \ c;
  beta = solve_tril(LS, R*z);
  invLsLp = solve_tril(LS, Lp');

  % Posterior mean and covariance
  f = invLsLp' * beta;
  Covf = invLsLp' * invLsLp;

end

% $$$ if issparse(Lp)
% $$$   density_in_gplearn = nnz(Lp) / prod(size(Lp));
% $$$ end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, df] = bound2cost(x, c, V, LV, inputs, kX, m, d, updatepseudo, ...
                              pseudoinputs)

if updatepseudo
  pseudoinputs = reshape(x(1:(m*d)), [d,m]);
  logtheta = x((m*d+1):end);
  [bound, dlogtheta, dpseudoinputs] = pseudobound(c, V, LV, inputs, ...
                                                  pseudoinputs, logtheta, kX);
else
  logtheta = x(:);
  [bound, dlogtheta] = pseudobound(c, V, LV, inputs, pseudoinputs, logtheta, ...
                                   kX);
end

f = -bound;

if updatepseudo
  df = -[dpseudoinputs(:); dlogtheta(:)];
else
  df = -dlogtheta(:);
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bound, dlogtheta, dpseudoinputs] = ...
    pseudobound(c, V, LV, inputs, pseudoinputs, logtheta, covfunc)

if ~isempty(pseudoinputs)
  %% PSEUDO INPUT CASE
  %disp('Pseudo case in gplearn > pseudobound');

  d = rows(pseudoinputs); % dimensionality of pseudo-inputs
  m = cols(pseudoinputs); % number of pseudo-inputs

  % GP-prior (co)variances and their gradients
  switch nargout
   case 1
    Kpp = feval(covfunc{:}, logtheta, pseudoinputs, pseudoinputs);
    Kxp = feval(covfunc{:}, logtheta, inputs, pseudoinputs);
    Kx = feval(covfunc{:}, logtheta, inputs, []);
   case 2
    [Kpp, dKpp_logtheta] = feval(covfunc{:}, logtheta, pseudoinputs, ...
                                   pseudoinputs);
    [Kxp, dKxp_logtheta] = feval(covfunc{:}, logtheta, inputs, pseudoinputs);
    [Kx, dKx_logtheta] = feval(covfunc{:}, logtheta, inputs, []);
   case 3
    [Kpp, dKpp_logtheta, dKpp_pseudo] = feval(covfunc{:}, logtheta, ...
                                                pseudoinputs, pseudoinputs);
    [Kxp, dKxp_logtheta, dKxp_pseudo] = feval(covfunc{:}, logtheta, ...
                                                inputs, pseudoinputs);
    [Kx, dKx_logtheta] = feval(covfunc{:}, logtheta, inputs, []);
   otherwise
    error('No defined action here')
  end

% $$$ if ~iscell(covfunc), covfunc = {covfunc}; end
% $$$ [K,tmp] = feval(covfunc{:}, logtheta, pseudoinputs, pseudoinputs);
% $$$ figure
% $$$ imagesc(K-Kpp)
% $$$ return

  Im = diag(sparse(ones(m,1)));
% $$$   norm_Kpp = normest(Kpp)
% $$$   Kpp = Kpp + 1e-14*norm_Kpp*Im; % SET PROPERLY
%  Kpp = Kpp + 1e-3*Im; % SET PROPERLY
% $$$   figure
% $$$   imagesc(Kpp);
% $$$   error('jou')
  
  % DEBUG: REGULARIZE
  %Kpp = regularize(Kpp);

  % Helpful variables
  Kpp = regularize(Kpp);
  try
    Lp = chol(Kpp, 'lower');
  catch
    warning('Covariance matrix ill conditioned');
    bound = -inf;
    dlogtheta = nan * zeros(length(logtheta),1);
    dpseudoinputs = nan * zeros(d,m);
    return
  end
%  [Lp,Kpp] = safechol(Kpp, 1e-10, 'lower');
  %Kpp = Lp*Lp';
  
  R = solve_tril(Lp, Kxp' * LV);
  S = Im + R*R';
  LS = chol(S, 'lower');
%  LS = safechol(S, 1e-16, 'lower');
  LA = Lp * LS;
  z = LV \ c;
  H = solve_tril(LS, R);
  T = solve_triu(LA', H);
  beta = H * z;
  b = T * z;
  invLpR = solve_triu(Lp', R);
  diagV = full(diag(V));
  dif = c - V*(Kxp*b);

  % Lower bound for loglikelihood
  bound = -logdettri(LS) - 0.5*z'*z + 0.5*beta'*beta ...
          - 0.5*Kx'*diagV ...
          + 0.5*traceprod(R,R,true);
  
  % General gradient coefficients
  Tpp = 0.5 * solve_triu(LA', solve_triu(LA', LS'*LS-eye(m))') ...
        - 0.5 * invLpR*invLpR';
  Txp = -LV * T' ...
        + LV*invLpR';

  % Gradient for hyperparameters
  if nargout >= 2
    n = length(logtheta);
    dlogtheta = zeros(n,1);
    for i=1:n
      dlogtheta(i) = traceprod(dKpp_logtheta(:,:,i), Tpp, true) ...
          - 0.5 * b' * dKpp_logtheta(:,:,i) * b ...
          + traceprod(dKxp_logtheta(:,:,i), Txp, true) ...
          + dif' * dKxp_logtheta(:,:,i) * b ...
          - 0.5 * dKx_logtheta(:,:,i)' * diagV;
    end
  end

  if nargout >= 3
    % Gradient for pseudo-inputs
    dpseudoinputs = zeros(d,m);
    Dpp = spalloc(m,m,2*m);
    for i=1:m
      for j=1:d
        Dpp(:) = 0;
        Dpp(i,:) = dKpp_pseudo(j,:,i);
        Dpp(:,i) = Dpp(:,i) + dKpp_pseudo(j,:,i)';
        
        dpseudoinputs(j,i) = traceprod(Dpp, Tpp, true) ...
            - 0.5 * b' * Dpp * b ...
            + dKxp_pseudo(j,:,i) * Txp(:,i) ...
            + dif' * dKxp_pseudo(j,:,i)' * b(i);
        
      end
    end
  end

  
else
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FULL CASE (NO PSEUDO INPUTS)
% $$$   disp('Full case in gplearn > pseudobound');
  d = rows(inputs); % dimensionality of pseudo-inputs
  n = cols(inputs);
  %m = cols(pseudoinputs); % number of pseudo-inputs

  % GP-prior (co)variances and their gradients
  switch nargout
   case 1
    K = feval(covfunc{:}, logtheta, inputs, inputs);
   case 2
    [K, dK_logtheta] = feval(covfunc{:}, logtheta, inputs, inputs);
   otherwise
    error('No pseudo inputs given');
  end

  In = diag(sparse(ones(n,1)));
% $$$   norm_K = normest(K);
% $$$   K = K + 1e-14*norm_K*In; % SET PROPERLY
%  K = K + 1e-6*In;
  
  % DEBUG: REGULARIZE
  %K = regularize(K);
  
  % ASSUME THAT V is DIAGONAL
  V = sparse(V);
  Ky = inv(V) + K;
  try
    L = chol(Ky, 'lower');
  catch
    warning('Covariance matrix ill conditioned');
    bound = -inf;
    dlogtheta = nan * zeros(length(logtheta),1);
    return
  end
%  [L,Ky] = safechol(Ky, 1e-10, 'lower');
  % DEBUG:
  %Ky = L*L';
  y = sparse(V) \ c;
  
  z = solve_tril(L, y);

  % Lower bound for loglikelihood
  bound = -logdettri(L) - 0.5*z'*z;
  
% $$$   % General gradient coefficients
% $$$   Tpp = 0.5 * solve_triu(LA', solve_triu(LA', LS'*LS-eye(m))') ...
% $$$         - 0.5 * invLpR*invLpR';
% $$$   Txp = -LV * T' ...
% $$$         + LV*invLpR';

  % Gradient for hyperparameters
  if nargout >= 2
    b = solve_triu(L',z);
    n = length(logtheta);
    dlogtheta = zeros(n,1);
    if issparse(Ky)
      invKy = sinv(Ky);
    else
      invKy = inv(Ky);
    end
    for i=1:n
      dlogtheta(i) = 0.5 * b' * dK_logtheta(:,:,i) * b ...
          - 0.5 * traceprod( invKy, dK_logtheta(:,:,i), true );
% $$$       dlogtheta(i) = traceprod(dKpp_logtheta(:,:,i), Tpp, true) ...
% $$$           - 0.5 * b' * dKpp_logtheta(:,:,i) * b ...
% $$$           + traceprod(dKxp_logtheta(:,:,i), Txp, true) ...
% $$$           + dif' * dKxp_logtheta(:,:,i) * b ...
% $$$           - 0.5 * dKx_logtheta(:,:,i)' * diagV;
    end
  end

end

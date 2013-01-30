% Q = vbfa(D, Y, W_module, X_module, noise_module, ...)
% 
% Variational Bayesian (VB) factor analysis (FA) learning algorithm with
% changeable modules for the latent variables.
%
% Optional parameters:
% 'maxiter'             100
% 'update_x'            1
% 'update_w'            1
% 'update_noise'        1
% 'rotate'              1
% 'rotation_checkgrad'  false
% 'rotation_show'       false
% 'rotation_maxiter'    10
% 'debug'               false
% 'autosave'            false
% 'autosavefile'        'vbfa_autosave'

function Q = vbfa(D, Y, W_module, X_module, noise_module, varargin)

[M,N] = size(Y);

options = struct('maxiter', 100, ...
                 'update_x', 1, ...
                 'update_w', 1, ...
                 'update_noise', 1, ...
                 'initialize_x', true, ...
                 'initialize_w', true, ...
                 'initialize_noise', true, ...
                 'rotate', 1, ...
                 'rotation_checkgrad', false, ...
                 'rotation_show', false, ...
                 'rotation_maxiter', 10, ...
                 'loglikelihood', true, ...
                 'debug', false, ...
                 'autosave', false,...
                 'autosavefile', 'vbfa_autosave');

[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Initialize
disp('Initializing variables..')
if options.initialize_w
  [Q.W,Q.CovW,Q.rhow] = W_module.initialize(D,M);
else
  QW = W_module.get_struct();
  Q.W = QW.X;
  Q.CovW = QW.CovX;
  Q.rhow = QW.rho;
end
if options.initialize_x
  [Q.X,Q.CovX,Q.rhox] = X_module.initialize(D,N);
else
  QX = W_module.get_struct();
  Q.X = QX.X;
  Q.CovX = QX.CovX;
  Q.rhox = QX.rho;
end
if options.initialize_noise
  [Q.Tau] = noise_module.initialize(M,N);
else
  QTau = noise_module.get_struct();
  Q.Tau = QTau.Tau;
end
%Q.W_module = W_module;
%Q.X_module = X_module;
%Q.noise_module = noise_module;


% Observed values
Obs = ~isnan(Y);
MN_obs = sum(Obs(:));

%%
%% VB learning

loglikelihood_old = -inf;
Q.loglikelihood = nan(options.maxiter,1);

Q.Y = Y;

% Replace missing values with zeros for computational reasons
Y(~Obs) = 0;

KL_W = -inf;
KL_X = -inf;
KL_Tau = -inf;

% Debugging stuff:
klw = nan(options.maxiter,1);
klx = nan(options.maxiter,1);
kltau = nan(options.maxiter,1);
loglike = nan(options.maxiter,1);


%f = zeros(4,options.maxiter);

disp('Starting the VB iteration..')
for ind=1:options.maxiter
  Q.ind = ind;
  
  startitercpu = cputime;
  
  %
  % Update variables
  %
  
  Tau = Q.Tau;
  Tau(~Obs) = 0;
  
  if index_selected(Q.ind, options.update_x)
    if index_selected(Q.ind, options.debug)
      disp('Update X')
    end
    %[U,Z] = projection(Y,Obs,W,CovW,Tau);
    [Q.X,Q.CovX,Q.rhox,logrhox,KL_X] = X_module.update(Q.ind, Y, Obs, Q.W, ...
                                                      Q.CovW, Q.rhow, Tau);
  end
  
  if index_selected(Q.ind, options.update_w)
    if index_selected(Q.ind, options.debug)
      disp('Update W')
    end
    %[U,Z] = projection(Y',Obs',X,CovX,Tau');
    [Q.W,Q.CovW,Q.rhow,logrhow,KL_W] = W_module.update(Q.ind, Y', Obs', Q.X, ...
                                                      Q.CovX, Q.rhox, Tau');
  end
  
  % Compute squared errors <(y_mn-w_m*x_n)^2>
  % (used by noise update and loglikelihood)
  E2 = zeros(size(Y));
  Yh = Q.W'*Q.X;
  E2 = Y.*Y;
  E2 = E2 - 2*Y.*Yh;
  E2 = spdiag(Q.rhow)*E2*spdiag(Q.rhox);
  if ndims(Q.CovW) == 2 && ndims(Q.CovX) == 2
    E2 = E2 ...
         + spdiag(Q.rhow)*(Yh.^2)*spdiag(Q.rhox) ...
         + spdiag(Q.rhow)*(Q.W'.^2)*Q.CovX ...
         + Q.CovW'*(Q.X.^2)*spdiag(Q.rhox) ...
         + Q.CovW'*Q.CovX;
  elseif ndims(Q.CovW) == 3 && ndims(Q.CovX) == 3
    xx = bsxfun(@times, reshape(Q.X*spdiag(Q.rhox),[D,1,N]), reshape(Q.X,[1,D,N])) ...
         + Q.CovX;
    xx = reshape(xx, [D*D,N]);
    ww = bsxfun(@times, reshape(Q.W*spdiag(Q.rhow),[D,1,M]), reshape(Q.W,[1,D,M])) ...
         + Q.CovW;
    ww = reshape(ww, [D*D,M]);
    E2 = E2 + ww'*xx;
  else
    % TODO: Optimize this..
    warning('Optimize this.. and is rho properly here??');
    for m=1:M
      for n=1:N
        if Obs(m,n)
          if ndims(Q.CovW) == 2
            ww = Q.W(:,m)*Q.W(:,m)' + diag(Q.CovW(:,m));
          else
            ww = Q.W(:,m)*Q.W(:,m)' + Q.CovW(:,:,m);
          end
          if ndims(Q.CovX) == 2
            xx = Q.X(:,n)*Q.X(:,n)' + diag(Q.CovX(:,n));
          else
            xx = Q.X(:,n)*Q.X(:,n)' + Q.CovX(:,:,n);
          end
          %WX_WX = WX_WX + traceprod(ww, xx);
          E2(m,n) = E2(m,n) + traceprod(ww, xx);
        end
      end
    end
  end
  E2(~Obs) = 0;

  if index_selected(Q.ind, options.update_noise)
    if index_selected(Q.ind, options.debug)
      disp('Update Tau')
    end
    
    [Q.Tau,LogTau,KL_Tau] = noise_module.update(Q.ind, E2, Obs);
%    [Tau,lowerbound] = noise_module.update(Y, Obs, v_W, W, CovW, v_X, X, CovX);
  end
  

  %
  % Rotate
  %
  
  if index_selected(Q.ind, options.rotate)
    
    % TODO: You could optimize the hyperparameters at the same time?
    disp('Rotating..')
    A = eye(D);
   
    if index_selected(Q.ind, options.rotation_checkgrad)
      mycheckgrad(@rotation_cost, A(:) + 0.5*randn(D^2,1), 1e-3, W_module, ...
                  X_module);
    end
    A = minimize(A(:), @rotation_cost, options.rotation_maxiter, W_module, ...
                 X_module);
    A = reshape(A, [D D]);
    if index_selected(Q.ind, options.rotation_show)
      A
    end
    [Q.W, Q.CovW] = W_module.rotate(A);
    [Q.X, Q.CovX] = X_module.rotate(inv(A)');

  end

  
  %
  % Evaluate VB lower bound
  %
  
  if index_selected(Q.ind, options.loglikelihood)
  % Likelihood part: <log p(Y|...)>
  logpdf_y = gaussian_logpdf(Q.Tau(Obs)'*E2(Obs), ...
                             0, ...
                             0, ...
                             -sum(LogTau(Obs))-logrhow'*sum(Obs,2)-sum(Obs,1)*logrhox, ...
                             MN_obs);

  % Debugging stuff
  klw(Q.ind) = -KL_W;
  klx(Q.ind) = -KL_X;
  kltau(Q.ind) = -KL_Tau;
  loglike(Q.ind) = logpdf_y;

  % Lower bound
  Q.loglikelihood(Q.ind) = logpdf_y - KL_W - KL_X - KL_Tau;
  
  if Q.loglikelihood(Q.ind) < loglikelihood_old
    warning(sprintf('Loglikelihood lower bound decreased relatively %e!', ...
                    (loglikelihood_old - Q.loglikelihood(Q.ind)) / ...
                    loglikelihood_old));
    % plot([Q.loglikelihood, klw, klx, kltau, loglike]);
  end
  
  fprintf('Iteration step %d: loglikelihood=%e (%.2f seconds)\n', Q.ind, ...
          Q.loglikelihood(Q.ind), cputime()-startitercpu);
  
  loglikelihood_old = Q.loglikelihood(Q.ind);
  end
  
  if index_selected(Q.ind, options.autosave)
    Q.W_struct = W_module.get_struct();
    Q.X_struct = X_module.get_struct();
    Q.Tau_struct = noise_module.get_struct();
    fprintf('Saving results to %s...', options.autosavefile);
    save(options.autosavefile, '-struct', 'Q');
    fprintf(' done.\n');
  end
  
end



function [c, dc] = rotation_cost(a, W_module, X_module)
N = sqrt(length(a));
A = reshape(a, N, N);
[U,S,V] = svd(A);
invS = diag(1./diag(S));
invAt = U*invS*V';
[c_w, dc_w] = W_module.rotation_cost(A, U, S, V);
[c_x, dc_x] = X_module.rotation_cost(invAt, U, invS, V);
dc_x = -invAt*dc_x'*invAt;
c = c_w + c_x;
dc = dc_w(:) + dc_x(:);


% Q = vbfa(D, Y, W_module, X_module, noise_module, ...)
% 
% Variational Bayesian (VB) factor analysis (FA) learning algorithm with
% changeable modules for the latent variables.

function Q = vbfa(D, Y, W_module, X_module, noise_module, varargin)

[M,N] = size(Y);

options = struct('maxiter', 100, ...
                 'update_x', 1, ...
                 'update_w', 1, ...
                 'update_noise', 1, ...
                 'rotate', 1, ...
                 'rotation_checkgrad', false, ...
                 'rotation_show', false, ...
                 'rotation_maxiter', 10, ...
                 'debug', false, ...
                 'autosave', false,...
                 'autosavefile', 'vbfa_autosave');

[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Initialize
disp('Initializing variables..')
[Q.v_W,Q.W,Q.CovW] = W_module.initialize();
[Q.v_X,Q.X,Q.CovX] = X_module.initialize();
[Q.Tau] = noise_module.initialize();
Q.W_module = W_module;
Q.X_module = X_module;
Q.noise_module = noise_module;


% Observed values
Obs = ~isnan(Y);
MN_obs = sum(Obs(:));

% Replace missing values with zeros for computational reasons
Y(~Obs) = 0;

%%
%% VB learning

loglikelihood_old = -inf;
Q.loglikelihood = nan(options.maxiter,1);

Q.Y = Y;

KL_W = -inf;
KL_X = -inf;
KL_Tau = -inf;

%f = zeros(4,options.maxiter);

disp('Starting the VB iteration..')
for ind=1:options.maxiter
  
  startitercpu = cputime;
  
  %
  % Update variables
  %
  
  if index_selected(ind, options.update_x)
    if index_selected(ind, options.debug)
      disp('Update X')
    end
    [Q.v_X,Q.X,Q.CovX,KL_X] = X_module.update(ind, Y, Obs, Q.v_W, Q.W, Q.CovW, Q.Tau);
  end
  
  if index_selected(ind, options.update_w)
    if options.debug, disp('Update W'), end
    if index_selected(ind, options.debug)
      disp('Update W')
    end
    [Q.v_W,Q.W,Q.CovW,KL_W] = W_module.update(ind, Y', Obs', Q.v_X, Q.X, Q.CovX, Q.Tau');
  end
  
  % Compute squared errors <(y_mn-w_m*x_n)^2>
  % (used by noise update and loglikelihood)
  E2 = zeros(size(Y));
  Yh = Q.W'*Q.X;
  E2 = Y.*Y;
  E2 = E2 - 2*Y.*Yh;
  if ndims(Q.CovW) == 2 && ndims(Q.CovX) == 2
    E2 = E2 + Yh.^2 + Q.W'.^2*Q.CovX + Q.CovW'*Q.X.^2 + Q.CovW'*Q.CovX;
  elseif ndims(Q.CovW) == 3 && ndims(Q.CovX) == 3
    xx = bsxfun(@times, reshape(Q.X,[D,1,N]), reshape(Q.X,[1,D,N])) + Q.CovX;
    xx = reshape(xx, [D*D,N]);
    ww = bsxfun(@times, reshape(Q.W,[D,1,M]), reshape(Q.W,[1,D,M])) + Q.CovW;
    ww = reshape(ww, [D*D,M]);
    E2 = E2 + ww'*xx;
  else
    % TODO: Optimize this..
    warning('Optimize this..');
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

  if index_selected(ind, options.update_noise)
    if index_selected(ind, options.debug)
      disp('Update Tau')
    end
    
    [Q.Tau,LogTau,KL_Tau] = noise_module.update(ind, E2, Obs);
%    [Tau,lowerbound] = noise_module.update(Y, Obs, v_W, W, CovW, v_X, X, CovX);
  end
  

  %
  % Rotate
  %
  
  if index_selected(ind, options.rotate)
    
    % TODO: You could optimize the hyperparameters at the same time?
    disp('Rotating..')
    A = eye(D);
   
    if index_selected(ind, options.rotation_checkgrad)
      mycheckgrad(@rotation_cost, A(:) + 0.5*randn(D^2,1), 1e-3, W_module, ...
                  X_module);
    end
    A = minimize(A(:), @rotation_cost, options.rotation_maxiter, W_module, ...
                 X_module);
    A = reshape(A, [D D]);
    if index_selected(ind, options.rotation_show)
      A
    end
    [Q.W, Q.CovW] = W_module.rotate(A);
    [Q.X, Q.CovX] = X_module.rotate(inv(A)');

  end

  
  %
  % Evaluate VB lower bound
  %
  
  % Likelihood part: <log p(Y|...)>
  logpdf_y = gaussian_logpdf(Q.Tau(Obs)'*E2(Obs), ...
                             0, ...
                             0, ...
                             -sum(LogTau(Obs)), ...
                             MN_obs);

  % Lower bound
  Q.loglikelihood(ind) = logpdf_y - KL_W - KL_X - KL_Tau;
  
% $$$   if ~isreal(Q.loglikelihood(ind))
% $$$     KL_W
% $$$     KL_X
% $$$     KL_Tau
% $$$     error('Loglikelihood imaginary')
% $$$   end
  
  if Q.loglikelihood(ind) < loglikelihood_old
    warning(sprintf('Loglikelihood lower bound decreased relatively %e!', ...
                    (loglikelihood_old - Q.loglikelihood(ind)) / loglikelihood_old));
  end
  
  fprintf('Iteration step %d: loglikelihood=%e (%.2f seconds)\n', ind, ...
          Q.loglikelihood(ind), cputime()-startitercpu);
  
% $$$   f(1,ind) = KL_W;
% $$$   f(2,ind) = KL_X;
% $$$   f(3,ind) = KL_Tau;
% $$$   f(4,ind) = loglikelihood(ind);
  
  loglikelihood_old = Q.loglikelihood(ind);
  
  if index_selected(ind, options.autosave)
    fprintf('Saving results to %s...', options.autosavefile);
    save(options.autosavefile, '-struct', 'Q');
    fprintf(' done.\n');
  end
  
end

% $$$ figure
% $$$ plot(f');

% $$$ if nargout >= 1
% $$$   
% $$$   Q.W_module = W_module;
% $$$   Q.W = W;
% $$$   Q.CovW = CovW;
% $$$   
% $$$   Q.X_module = X_module;
% $$$   Q.X = X;
% $$$   Q.CovX = CovX;
% $$$   
% $$$   Q.noise_module = noise_module;
% $$$   Q.Tau = Tau;
% $$$   
% $$$   Q.loglikelihood = loglikelihood;
% $$$   
% $$$ end

%function rotate(W_module, X_module)
%end

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

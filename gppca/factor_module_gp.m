% GP prior for the components. 
%
% module = factor_module_gp(N, covfuncs, theta, varargin)

function module = factor_module_gp(N, covfuncs, theta, varargin)

error('Not yet implemented')

D = length(covfuncs);

options = struct('is_pseudo', false(D,1), ...
                 'noise_prior', [], ...
                 'update_hyperparameters', true, ...
                 'maxiter_hyperparameters', 30, ...
                 'init', []);

[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

K_xx = cell(D,1);
K_px = cell(D,1);
K_pp = cell(D,1);

rho = [];
X = [];
CovX = [];

linsolve_K = cell(D,1);

module.get_struct = @get_struct;
module.initialize = @initialize;
module.update = @update;
module.rotate = @rotate;
module.rotation_cost = @rotation_cost;

  function S = get_struct()
  % TODO: not ready yet..
  S.X = X;
  S.CovX = CovX;
  S.rho = rho;
  S.options = options;
  S.covfuncs = covfuncs;
  S.theta = theta;
  end

  function [Z, CovZ, rhoz] = initialize()
  if isempty(options.init)
    X = zeros(D,N);
    for d=1:D
      if ~options.is_pseudo(d)
        K_xx{d} = covfuncs{d}(theta{d});
        L = chol(K_xx{d}, 'lower');
        X(d,:) = L * randn(N,1);
      else
        [K_pp{d}, K_px{d}, K_xx{d}] = covfuncs{d}(theta{d});
        L = chol(K_pp{d}, 'lower');
        Xp = L * randn(length(L),1);
        X(d,:) = K_px{d}' * linsolve_lchol(L, Xp);
      end
    end
  else
    X = options.init;
  end
  CovX = 1*ones(D,N);
  rho = ones(N,1);
  
  Z = X;
  CovZ = CovX;
  rhoz = rho;
  end

  function [Z, CovZ, rhoz, logrhoz, KL] = update(iter, Y, Obs, W, CovW, rhow, Tau)

  linsolve_K(:) = {[]};
  
  M = size(Tau,1);
  
  % Update each component (component-wise factorization)
  KL = 0;
  Tau(~Obs) = 0;
  for d=1:D
    
    % Indeces of the other components
    others = [1:(d-1), (d+1):D];

    %
    % Pre-evaluations
    %
    
    % TODO: What if there are totally missing rows in Y??
    % Then u will contain zeros..
    
    if ndims(CovW) == 2 && all(size(CovW)==[D,M])
      Yh = W(others,:)'*X(others,:);
      u = Tau' * (W(d,:).^2+CovW(d,:))';
      z = (Tau.*(Y-Yh))' * W(d,:)';
    else
      % TODO: Optimize this.. Do not use loops..
      %
      % See help from ARD-module..
      warning('Optimize this..');
      u = zeros(N,1);
      z = zeros(N,1);
      for n=1:N
        obs = Obs(:,n);
          
        % Is this correct
        covw = reshape(CovW(d,others,obs),[length(others), sum(obs)]);
        z(n) = (Tau(obs,n).*v_W(obs).*W(d,obs)')' * (Y(obs,n)-W(others,obs)) - ...
               (Tau(obs,n)'*covw')*X(others,n);

        varw = vec(CovW(d,d,obs));
        ww = v_W(obs).*(W(d,obs)'.^2) + varw;
        u(n) = vec(Tau(obs,n))' * ww;
      end
    end
    % Remove totally missing observations
    ind_obs = (u~=0);
    z = z(ind_obs);
    u = u(ind_obs);
    % Helpful variables
    invU_z = z./u;
    invU = spdiag(1./u);
    % Noise covariance function
    covfunc_noise = gp_cov_wrap(invU);
    
    %
    % Update variables
    %
    
    if ~options.is_pseudo(d)
      
      % Standard full GP
      
      % Form the joint covariance function
      covfunc = gp_cov_sum(gp_cov_select(covfuncs{d}, ind_obs), covfunc_noise);
      
      % Learn the hyperparameters
      if index_selected(iter, options.update_hyperparameters)
        [theta{d}, logpdf] = gp_learn(invU_z, covfunc, theta{d}, ...
                                      'maxiter', options.maxiter_hyperparameters, ...
                                      'checkgrad', false);
        K_xx{d} = covfuncs{d}(theta{d});
      else
        if isempty(K_xx{d})
          K_xx{d} = covfuncs{d}(theta{d});
        end
        logpdf = gp_loglikelihood(invU_z, covfunc, theta{d});
      end
      
      % Evaluate the mean and the variance
      y = invU_z;
      K_y = K_xx{d}(ind_obs,ind_obs) + invU;
      K_yf = K_xx{d}(ind_obs,:);
      k_f = diag(K_xx{d});
      [X(d,:),CovX(d,:)] = gp_predict(y, K_y, K_yf, k_f);
% $$$       [X(d,:),CovX(d,:)] = gp_predict(invU_z, K_xx{d}+invU, K_xx{d}, ...
% $$$                                  diag(K_xx{d}));
      
    else
      
      % Sparse approximation GP using pseudo inputs
      
      % Remove totally missing observations
      covfunc = gp_cov_select_pseudo(covfuncs{d}, ind_obs);
      
      % Learn the hyperparameters
      if index_selected(iter, options.update_hyperparameters)

        [theta{d}, ~, logpdf] = gp_learn_pseudo(invU_z, ...
                                                covfunc, ...
                                                theta{d}, ...
                                                covfunc_noise, ...
                                                [], ...
                                                'maxiter', ...
                                                options.maxiter_hyperparameters, ...
                                                'checkgrad', false);
        [K_pp{d}, K_px{d}, K_xx{d}] = covfuncs{d}(theta{d});
      else
        logpdf = gp_loglikelihood_pseudo(invU_z, covfunc, theta{d}, ...
                                         covfunc_noise, []);
        if isempty(K_pp{d}) || isempty(K_px{d}) || isempty(K_xx{d})
          [K_pp{d}, K_px{d}, K_xx{d}] = covfuncs{d}(theta{d});
        end
      end
      
      % Evaluate the mean and the variance
      [X(d,:),CovX(d,:)] = gp_predict_pseudo(invU_z, ...
                                             @(A) bsxfun(@times, u, A), ...
                                             K_px{d}(:,ind_obs), ...
                                             K_pp{d}, ...
                                             K_px{d}, ...
                                             K_xx{d});
      
    end
    
    %
    % Compute Kullback-Leibler divergence KL(q(X)||p(X))
    %
    logpdf_bound = gaussian_logpdf(z'*invU_z, ...
                                   z'*X(d,ind_obs)', ...
                                   (X(d,ind_obs).^2+CovX(d,ind_obs))*u, ...
                                   -sum(log(u)), ...
                                   N);
    
    KL = KL + logpdf_bound - logpdf;
    
  end
  
  % Update q(rho)
  rho = ones(N,1);
  logrho = zeros(N,1);
  KL_rho = 0;
  
  KL = KL + KL_rho;
  
  Z = X;
  CovZ = CovX;
  rhoz = rho;
  logrhoz = logrho;

  end
  
  function [c, dc] = rotation_cost(A, U, S, V)
% $$$   if any(options.is_pseudo)
% $$$     error('Rotation cost not implemented for pseudo inputs yet..');
% $$$   end
  % Cost from <log q(X)>
  c = -N * logdet_diag(S);
  dc = -N * inv(A)'; % FIX THIS
  % Cost from -<log p(X)>
  AX = A*X;
  varAX = A.^2*CovX;
  invdiagK = zeros(D,N);
  invK_XA = zeros(N,D);
  for d = 1:D
    if ~options.is_pseudo(d)
      if isempty(linsolve_K{d})
        linsolve_K{d} = get_linsolve_cov(K_xx{d});
      end
      invdiagK(d,:) = 1./diag(K_xx{d});
      invK_XA(:,d) = linsolve_K{d}(AX(d,:)');
    else
      %error('Pseudo cost not maybe correct..');
      if isempty(linsolve_K{d})
        L = chol(K_pp{d}, 'lower');
        Z = linsolve_tril(L, K_px{d});
        h = K_xx{d} - dot(Z,Z,1)';
        h = h + 1e-6*max(h);
        %H = spdiag(h);
        %Z = L(K_xp{d}*diag(h.^(-0.5)));
        H = spdiag(h);
        %plot(h)
        Lambda = speye(size(L)) + Z*spdiag(1./h)*Z';
        L_Lambda = get_linsolve_cov(Lambda);
        linsolve_K{d} = @(x) (H\x - H\(Z'*L_Lambda(Z*(H\x))));
        
% $$$         L = get_linsolve_cov(K_pp{d});
% $$$         h = K_xx{d} - diagprod(K_xp{d}, L(K_xp{d}')');
        %H = 1e-6*speye(size(H)); % better conditioning...
        %LL = get_linsolve_cov(K_pp{d}+K_xp{d}'*(H\K_xp{d}));
        %linsolve_K{d} = @(x) (H\x - H\(K_xp{d}*LL(K_xp{d}')*(H\x)));
      end
      invK_XA(:,d) = linsolve_K{d}(AX(d,:)');
      invdiagK(d,:) = 1./K_xx{d};
      %invK_XA(:,d) = AX(d,:)' ./ K_xx{d};
    end
    c = c + 0.5 * v_X * (AX(d,:)*invK_XA(:,d));
    c = c + 0.5 * varAX(d,:) * invdiagK(d,:)';
  end
  dc = dc + v_X * (invK_XA'*X');
  dc = dc + A .* (invdiagK*CovX');
% $$$   dc = dc + 0.5 * (invdiagK*CovX'*A' + A'*invdiagK*CovX');
  
  end
  
  function [Z, CovZ] = rotate(A)
  X = A*X;
  CovX = A.^2*CovX;
  Z = X;
  CovZ = CovX;
  end

end

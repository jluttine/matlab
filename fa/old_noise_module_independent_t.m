% module = noise_module_independent_t(noise_module, varargin)
%
% nu can be 'pooled', 'separate_rows', 'separate_columns'

function module = noise_module_independent_t(varargin)


options = struct('init', struct(), ...
                 'nu', 'pooled', ...
                 'update_nu', true, ...
                 'update_u', true);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

a_U = [];
b_U = [];
U = [];
nu = [];

module.initialize = @initialize;
module.update = @update;
module.get_struct = @get_struct;

  function [Tau] = initialize(M, N)

  % Allocate memory
  switch options.nu
   case 'pooled'
    nu = NaN;
   case 'separate_rows'
    nu = nan(M,1);
   case 'separate_columns'
    nu = nan(1,N);
  end
  U = nan(M,N);
  a_U = nan(M,N);
  b_U = nan(M,N);

  % Parse custom initialization
  init = struct('U', 1, ...
                'nu', 1);
  [init, errmsg] = argparse(init, options.init);
  error(errmsg);
  
  % Initialize
  U(:,:) = init.U;
  nu(:) = init.nu;
  Tau = U;
  end
  
  function S = get_struct()
  % TODO: not ready yet..
  S.module = 'noise_module_independent_t';
  S.nu = nu;
  S.U = U;
  S.a_U = a_U;
  S.b_U = b_U;
  S.options = options;
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  % OBS is the number of observations for each element of E2.
  
  [M,N] = size(E2);
  
  %
  % Update nu
  %
  
  if index_selected(iter, options.update_nu)
    switch options.nu
     
     case 'pooled'
      nu = t_ml(E2(Obs), nu(1));
     
     case 'separate_rows'
      for m=1:M
        nu(m) = t_ml(E2(m,Obs(m,:)), nu(m));
      end
      %lognu_debug = log(nu)
     
     case 'separate_columns'
      for n=1:N
        nu(n) = t_ml(E2(Obs(:,n),n), nu(n));
      end
     
     otherwise
      error('Unknown model for nu');
      
    end
  end
  
  %
  % Update U
  %
  
  if index_selected(iter, options.update_u)
    
    % Update the distribution
    a_U = bsxfun(@plus, nu, Obs) / 2;
    b_U = bsxfun(@plus, nu, E2) / 2;
    U = a_U ./ b_U;
    LogU = psi(a_U) - log(b_U);
    
    % KL-term: <log q(U)> - <log p(U)>
    KL_U = 0;
    KL_U = KL_U - sum(gamma_entropy(a_U(Obs),b_U(Obs)));
    switch options.nu
      
     case 'pooled'
      KL_U = KL_U - sum(gamma_logpdf(U(Obs), nu/2, nu/2, LogU(Obs)));
      
     case 'separate_rows'
      for m=1:M
        Nmv = Obs(m,:);
        KL_U = KL_U - sum(gamma_logpdf(U(m,Nmv),nu(m)/2,nu(m)/2,LogU(m,Nmv)));
      end
      
     case 'separate_columns'
      for n=1:N
        Mmv = Obs(:,n);
        KL_U = KL_U - sum(gamma_logpdf(U(Mmv,n),nu(n)/2,nu(n)/2,LogU(Mmv,n)));
      end
      
     otherwise
      error('Unknown model for nu');
      
    end
    
  else
    
    KL_U = NaN;
    
  end    

  Tau = U;
  LogTau = LogU;
  KL_Tau = KL_U;
  
  
  end

end

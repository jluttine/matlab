% NOISE_MODULE_INDEPENDENT_T  -  Independent Student-t noise module for
%                                variational Bayesian factor models
%
% MODULE = NOISE_MODULE_INDEPENDENT_T(M, N, ...)
%
% This module models only the robustness of each sample but not the noise
% level.  Use in combination with other non-robust noise modules.
%
% M : Rows in the data matrix
% N : Columns in the data matrix
%
% Optional parameters:
% 'nu'   : Tells whether the degrees of freedom is common for all
%          elements or may vary between rows/columns:
%            'pooled' : One common degrees of freedom parameter
%            'separate_rows' : Separate degrees of freedom parameter for
%                              each row of the data matrix
%            'separate_columns' : Separate degrees of freedom parameter
%                                 for each column of the data matrix
% 'init' : A struct containing initial values for the degrees of
%          freedom and the posterior Gamma distributions of U(J):
%            'nu'  : Degrees of freedom; scalar or Mx1 vector or 1xN
%                    vector (default: 1) 
%            'a_U' : Shape parameters; MxN matrix (default: 0.5)
%            'b_U' : Scale parameter; MxN matrix (default: 0.5)

% Last modified 2011-08-22
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)


function module = noise_module_independent_t(M,N,varargin)


options = struct('init', struct(), ...
                 'nu', 'pooled', ...
                 'update_nu', true, ...
                 'update_u', true);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Allocate memory
switch options.nu
 case 'pooled'
  nu = NaN;
 case 'separate_rows'
  nu = nan(M,1);
 case 'separate_columns'
  nu = nan(1,N);
end
a_U = nan(M,N);
b_U = nan(M,N);

% Parse custom initialization
init = struct('a_U', 1/2, ...
              'b_U', 1/2, ...
              'nu', 1);
[init, errmsg] = argparse(init, options.init);
error(errmsg);

% Initialize
nu(:) = init.nu;
a_U(:,:) = init.a_U;
b_U(:,:) = init.b_U;
U = a_U ./ b_U;
LogU = psi(a_U) - log(b_U);

module.update = @update;
module.get_struct = @get_struct;

  function S = get_struct()
  S.module = 'noise_module_independent_t';
  S.Tau = U;
  S.posterior = struct('nu', nu, ...
                       'a_U', a_U, ...
                       'b_U', b_U);
  S.init = init;
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

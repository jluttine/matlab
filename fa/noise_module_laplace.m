% NOISE_MODULE_LAPLACE  -  Laplacian noise module for variational
%                          Bayesian factor models
%
% MODULE = NOISE_MODULE_LAPLACE(M, N, ...)
%
% M : Rows in the data matrix
% N : Columns in the data matrix
%
% Optional parameters:
% 'init'  : A struct containing initial values:
%             'U' : MxN matrix of initial precisions (default: 1)

% Last modified 2011-08-16
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@aalto.fi)

function module = noise_module_laplace(M,N,varargin)


options = struct('init', struct(), ...
                 'update_u', true);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Parse custom initialization
init = struct('U', 1);
[init, errmsg] = argparse(init, options.init);
error(errmsg);

% Initialize
U = nan(M,N);
LogU = nan(M,N);
U(:,:) = init.U;

module.update = @update;
module.get_struct = @get_struct;

  function S = get_struct()
  S.module = 'noise_module_independent_t';
  S.Tau = U;
  S.posterior = struct('U', U);
  S.init = init;
  S.options = options;
  end

  function [Tau, LogTau, KL_Tau] = update(iter, E2, Obs)

  % OBS is the number of observations for each element of E2.
  
  [M,N] = size(E2);
  
  %
  % Update U
  %
  
  if index_selected(iter, options.update_u)
    
    % Update the distribution
    
    % a_U = bsxfun(@plus, nu, Obs) / 2;
    % b_U = bsxfun(@plus, nu, E2) / 2;
    % U = a_U ./ b_U;
    U(Obs) = E2(Obs) .^ (-0.5);
    LogU = nan(M,N); %psi(a_U) - log(b_U);
  end
  
% $$$   if iter == 4
% $$$     figure
% $$$     imagesc(U);
% $$$     error('debuggin')
% $$$   end
% $$$   if iter == 1
% $$$     figure
% $$$     imagesc(U);
% $$$     error('debuggin')
% $$$   end
  
  %minU = min(U(:));
  %maxU = max(U(:));
  
  % KL-term: <log q(U)> - <log p(U)>
  KL_U = NaN;
% $$$     0;
% $$$     KL_U = KL_U - sum(gamma_entropy(a_U(Obs),b_U(Obs)));
% $$$     switch options.nu
% $$$       
% $$$      case 'pooled'
% $$$       KL_U = KL_U - sum(gamma_logpdf(U(Obs), nu/2, nu/2, LogU(Obs)));
% $$$       
% $$$      case 'separate_rows'
% $$$       for m=1:M
% $$$         Nmv = Obs(m,:);
% $$$         KL_U = KL_U - sum(gamma_logpdf(U(m,Nmv),LogU(m,Nmv),nu(m)/2,nu(m)/2));
% $$$       end
% $$$       
% $$$      case 'separate_columns'
% $$$       error('Not implemented yet')
% $$$       
% $$$      otherwise
% $$$       error('Unknown model for nu');
% $$$       
% $$$     end
% $$$     
% $$$   else
% $$$     
% $$$     KL_U = NaN;
% $$$     
% $$$   end    

  %
  % Compute the joint noise stuff
  %
  
  Tau = U;
  LogTau = LogU;
  KL_Tau = KL_U;
  
  
  end

end

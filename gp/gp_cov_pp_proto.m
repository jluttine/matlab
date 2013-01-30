
% GP_COV_PP - Piecewise polynomial covariance function with compact support
% for Gaussian processes.
%
% [K, DK_LOGTHETA, DK_X2] = GP_COV_CS(X1, X2, LOGTHETA)

% Last modified 2010-10-28
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_pp_proto(x1,x2)

if nargin < 2
  x2 = x1;
end

% Dimensions
d = rows(x2); % dimensionality of inputs
m = cols(x1); % number of other inputs
n = cols(x2); % number of inputs

q = 2; % PP parameter, fix to 2 for now..

% Distance matrix
D = sqrt(sq_dist(x1,x2));

covfunc = @get_covariance;

  function [K,L,M] = get_covariance(theta)
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    K = 1;
    L = m;
    M = n;
    return
  end
  
  % Threshold
  thres = theta(1); % length scale, threshold
  R = D./thres;
  I = sparse(R<1);
  NNZ = nnz(I); % number of nonzero elements

  j = floor(rows(x1)) + q + 1;
  switch q
    %%%%%%%%%%
   case 0
    error('not yet implemented for q=0')
    K(I) = (1-R(I)) .^ j;
    %%%%%%%%%%
   case 1
    error('not yet implemented for q=1')
    %%%%%%%%%%
    
   case 2
    a = j^2 + 4*j + 3;
    b = 3*j + 6;
    c = 3;
    z = 1/3;
    
    % Helper
    R2 = spalloc(m,n,NNZ);
    R2(I) = D(I).^2 ./ (thres^2);

    % Covariance matrix
    K = spalloc(m,n,NNZ);
    K(I) = z * (1-R(I)).^(j+2) .* (a*R2(I) + b*R(I) + c);
    
    dK_dR = spalloc(m,n,NNZ); % copy sparseness
    dK_dR(I) = -z * (j+2) * (1-R(I)).^(j+1) .* (a*R2(I) + b*R(I) + c) ...
        + z * (1-R(I)).^(j+2) .* (2*a*R(I) + b);

    % Gradient for hyperparameters
    if nargout >= 2
      % TODO 3D SPARSE MATRICES DO NOT WORK. USE CELL ARRAYS?
      dK_dlogtheta = spalloc(m,n,NNZ);
      dK_dlogtheta(I) = dK_dR(I) .* (-D(I) ./ thres);
      dK_dlogtheta = full(dK_dlogtheta); % blaaah.. :(
    end

    % Gradients for inputs x2
    if nargout >= 3
      if isempty(x2)
        error('Can''t calculate gradient: x2 not given');
      end
      % TODO 3D SPARSE MATRICES DO NOT WORK. USE CELL ARRAYS?
      % TODO: you have dK_dR and dD_dx but NOT dR_dD!!
      dK_dx2 = bsxfun(@times, reshape(full(dK_dR./thres),[1,m,n]), dD_dx2);
      % blaah the need for fullness.. :(
    end
    
    %%%%%%%%%%
    
   case 3
    error('not yet implemented for q=3')
    %%%%%%%%%%
   otherwise
    error('q not valid')
  end


  if ~issparse(K)
    error('K not sparse, it should, wtf?!');
  end
  
  if nargout >= 2
    L = ldlchol2chol( ldlchol(K) );
  end
  end
  
end


% GP_COV_PP - Piecewise polynomial covariance function with compact support
%             for Gaussian processes.
%
%   FUNC = GP_COV_PP(D, DIM)
%
% D is a distance matrix and DIM is the dimensionality of the (Euclidean)
% input space.
%
%   K = FUNC(THETA)
%   [K, DK] = FUNC(THETA)
%
% THETA is the cut-off distance.

% Last modified 2010-12-02
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function covfunc = gp_cov_pp(D, d)

error(nargchk(2,inf,nargin,'struct'));

% PP parameter, fix to 2 for now..
q = 2;

covfunc = @get_covariance;

  function varargout = get_covariance(theta)
  
  varargout = cell(nargout,1);
  
  if nargin == 0
    % Return the number of parameters and the size of the covariance matrix
    if nargout >= 1
      varargout{1} = 1; % number of parameters
      if nargout >= 2
        varargout{2} = size(D,1); % number of rows
        if nargout == 3
          varargout{3} = size(D,2); % number of columns
        else
          error('Too many outputs');
        end
      end
    end
    return
  end
  
  if numel(theta) ~= 1
    error('Wrong number of parameters (%d), should be 1', numel(theta));
  end
  
  % Threshold
  thres = theta(1); % length scale, threshold
  R = D./thres;
  I = sparse(R<1);
  NNZ = nnz(I); % number of nonzero elements

  j = floor(d/2) + q + 1;
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
    R2 = spalloc(size(D,1),size(D,2),NNZ);
    R2(I) = R(I).^2; % D(I).^2 ./ (thres^2);

    % Covariance matrix
    K = spalloc(size(D,1),size(D,2),NNZ);
    K(I) = z * (1-R(I)).^(j+2) .* (a*R2(I) + b*R(I) + c);
    varargout{1} = K;
    
    % Gradient for hyperparameters
    if nargout >= 2
      dK_dR = spalloc(size(D,1),size(D,2),NNZ); % copy sparseness
      dK_dR(I) = -z * (j+2) * (1-R(I)).^(j+1) .* (a*R2(I) + b*R(I) + c) ...
          + z * (1-R(I)).^(j+2) .* (2*a*R(I) + b);

      dK_dlogtheta = spalloc(size(D,1),size(D,2),NNZ);
      dK_dlogtheta(I) = dK_dR(I) .* (-D(I) .* thres^(-2));
      varargout{2} = {dK_dlogtheta};
      % TODO 3D SPARSE MATRICES DO NOT WORK. USE CELL ARRAYS?
% $$$       dK_dlogtheta = full(dK_dlogtheta); % blaaah.. :(
    end

% $$$     % Gradients for inputs x2
% $$$     if nargout >= 3
% $$$       if isempty(x2)
% $$$         error('Can''t calculate gradient: x2 not given');
% $$$       end
% $$$       % TODO 3D SPARSE MATRICES DO NOT WORK. USE CELL ARRAYS?
% $$$       % TODO: you have dK_dR and dD_dx but NOT dR_dD!!
% $$$       dK_dx2 = bsxfun(@times, reshape(full(dK_dR./thres),[1,m,n]), dD_dx2);
% $$$       % blaah the need for fullness.. :(
% $$$     end
    
    %%%%%%%%%%
    
   case 3
    error('not yet implemented for q=3')
    %%%%%%%%%%
   otherwise
    error('q not valid')
  end

  end
  
end


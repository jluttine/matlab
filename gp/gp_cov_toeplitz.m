% Uses SPTOEPLITZ for sparse matrices
function covfunc_toeplitz = gp_cov_toeplitz(covfunc)
covfunc_toeplitz = @cov_toeplitz;
  function varargout = cov_toeplitz(theta)
  varargout = cell(nargout,1);
  if nargin == 0
    if nargout == 1
      varargout{1} = covfunc();
    else
      out = cell(max(nargout,3),1);
      [out{:}] = covfunc();
      out{2} = max(out{2}, out{3});
      out{3} = out{2};
      varargout(:) = out(1:nargout);
    end
    return
  end
  if nargout == 1
    k = covfunc(theta);
    if issparse(k)
      varargout{1} = sptoeplitz(k);
    else
      varargout{1} = toeplitz(k);
    end
  elseif nargout == 2
    [k, dk] = covfunc(theta);
    if issparse(k)
      varargout{1} = sptoeplitz(k);
    else
      varargout{1} = toeplitz(k);
    end
    varargout{2} = cell(numel(dk),1);
    for i=1:numel(dk)
      if issparse(k)
        varargout{2}{i} = sptoeplitz(dk{i});
      else
        varargout{2}{i} = toeplitz(dk{i});
      end
    end
  end
  end
end


% Uses SPTOEPLITZ for sparse matrices
function covfunc_toeplitz = gp_cov_toeplitz_block(covfunc)
covfunc_toeplitz = @cov_toeplitz_block;
  function varargout = cov_toeplitz_block(theta)
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
    varargout{1} = toeplitz_block(k);
  elseif nargout == 2
    [k, dk] = covfunc(theta);
    varargout{1} = toeplitz_block(k);
    varargout{2} = cell(numel(dk),1);
    for i=1:numel(dk)
      varargout{2}{i} = toeplitz_block(dk{i});
    end
  end
  end
end


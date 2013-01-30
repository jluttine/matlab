function func = gp_cov_select(covfunc, ind1, ind2)
if nargin < 3
  ind2 = ind1;
end
func = @cov_select;
  function varargout = cov_select(theta)
  out = cell(nargout,1);
  varargout = cell(nargout,1);
  if nargin == 0
    [out{:}] = covfunc();
    varargout{1} = out{1};
    if islogical(ind1)
      varargout{2} = sum(ind1);
    elseif isvector(ind1)
      varargout{2} = numel(ind1);
    else
      error('ind1 should be a vector')
    end
    if islogical(ind2)
      varargout{3} = sum(ind2);
    elseif isvector(ind2)
      varargout{3} = numel(ind2);
    else
      error('ind2 should be a vector')
    end
    varargout(4:end) = out(4:end);
    return
  end
  [out{:}] = covfunc(theta);
  varargout{1} = out{1}(ind1,ind2);
  if nargout >= 2
    varargout{2} = cell(numel(theta),1);
    for n=1:numel(theta)
      % Scale the gradients
      varargout{2}{n} = out{2}{n}(ind1,ind2);
    end
  end
  end
end

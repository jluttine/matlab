function func = gp_cov_remove(covfunc, ind1, ind2)
if nargin < 3
  ind2 = ind1;
end
func = @cov_remove;
  function varargout = cov_remove(theta)
  ind1
  out = cell(nargout,1);
  varargout = cell(nargout,1);
  if nargin == 0
    [out{:}] = covfunc();
    varargout{1} = out{1};
    varargout(2:end) = out(2:end);
    return
  end
  [out{:}] = covfunc(theta);
  out{1}(ind1,ind2) = [];
  varargout{1} = out{1};
  if nargout >= 2
    varargout{2} = cell(numel(theta),1);
    for n=1:numel(theta)
      % Scale the gradients
      out{2}{n}(ind1,ind2) = [];
      varargout{2}{n} = out{2}{n};
    end
  end
  end
end

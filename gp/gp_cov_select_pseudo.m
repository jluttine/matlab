function func = gp_cov_select_pseudo(covfunc, ind)
func = @cov_select;
  function varargout = cov_select(theta)
  out = cell(nargout,1);
  varargout = cell(nargout,1);
  if nargin == 0
    error('This feature not yet implemented');
% $$$     [out{:}] = covfunc();
% $$$     varargout{1} = out{1};
% $$$     if islogical(ind1)
% $$$       varargout{2} = sum(ind1);
% $$$     elseif isvector(ind1)
% $$$       varargout{2} = numel(ind1);
% $$$     else
% $$$       error('ind1 should be a vector')
% $$$     end
% $$$     if islogical(ind2)
% $$$       varargout{3} = sum(ind2);
% $$$     elseif isvector(ind2)
% $$$       varargout{3} = numel(ind2);
% $$$     else
% $$$       error('ind2 should be a vector')
% $$$     end
% $$$     varargout(4:end) = out(4:end);
    return
  end
  [out{:}] = covfunc(theta);
  if nargout >= 1
    varargout{1} = out{1};
  end
  if nargout >= 2
    varargout{2} = out{2}(:,ind);
  end
  if nargout >= 3
    varargout{3} = out{3}(ind);
  end
  if nargout >= 4
    varargout{4} = out{4};
  end
  if nargout >= 5
    varargout{5} = cell(numel(theta),1);
    for n=1:numel(theta)
      % Scale the gradients
      varargout{5}{n} = out{5}{n}(:,ind);
    end
  end
  if nargout >= 6
    varargout{6} = cell(numel(theta),1);
    for n=1:numel(theta)
      % Scale the gradients
      varargout{6}{n} = out{6}{n}(ind);
    end
  end
  end
end

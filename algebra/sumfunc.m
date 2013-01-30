function g = sumfunc(f1, f2)
g = @summing;
  function varargout = summing(varargin)
  nout = max(nargout,1);
  out1 = cell(1,nout);
  [out1{:}] = feval(f1,varargin{:});
  out2 = cell(1,nout);
  [out2{:}] = feval(f2,varargin{:});
  varargout = cellsum(out1,out2);
  end
end

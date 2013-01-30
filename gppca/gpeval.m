function gpeval(func, varargin)


if isa(func, 'function_handler');
  f = func2str(func);
end
  
eval(sprintf('feval(%s, varargin{:})', f));

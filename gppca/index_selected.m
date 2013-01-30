function bool = index_selected(iter, indeces)

if islogical(indeces)
  if isscalar(indeces)
    bool = indeces;
  else
    error('Logical values must not be matrices');
  end
elseif isscalar(indeces)
  if indeces > 0
    bool = iter >= indeces;
  else
    error('Index must be positive integer.');
  end
else
  bool = any(iter==indeces);
end

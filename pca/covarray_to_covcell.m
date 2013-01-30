function Covcell = covarray_to_covcell(Covarray)

if isnumeric(Covarray)
  [M,N,D] = size(Covarray);
  Covcell = cell(D,1);
  for d=1:D
    Covcell{d} = Covarray(:,:,d);
  end
end
function Covarray = covcell_to_covarray(Covcell)

if iscell(Covcell)
  sz = size(Covcell{1});
  N = length(Covcell);
  Covarray = zeros([sz, N]);
  for n=1:N
    Covarray(:,:,n) = Covcell{n};
  end
end

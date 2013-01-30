function func = get_linsolve_cov(K)
if issparse(K)
  [LD,p,q] = ldlchol(K);
  if p ~= 0
    figure
    imagesc(K);
    error('Matrix not positive definite')
  end
  func = @(x) linsolve_ldlchol(LD, x, q);
else
  [L,p] = chol(K, 'lower');
  if p ~= 0
    figure
    imagesc(K);
    error('Matrix not positive definite')
  end
  func = @(x) linsolve_lchol(L, x);
end


function x = mymvnrnd(mu, Cov)
warning('This function is deprecated')


if false %issparse(Cov)
  [L,p,S] = chol(Cov, 'lower');
  x = S' * L * randn(rows(Cov),1) + mu;
else
  L = chol(Cov, 'lower');
  x = L * randn(rows(Cov),1) + mu;
end
function y = gprnd(X, logtheta, covfunc)

if ~iscell(covfunc)
  covfunc = {covfunc};
end

K = regularize(feval(covfunc{:}, logtheta, X, X));% + 1e-5*diag(sparse(ones(cols(X),1)));

y = mymvnrnd(0, K);
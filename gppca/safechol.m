function [U, Xreg] = safechol(X, coeff, varargin)
% function [U,X] = safechol(X,coeff,...)
%
% Returns the Cholesky decomposition such that U'*U = X.
% X must be symmetric positive definite matrix. Tries to solve numerical
% issues by adding a small value to the diagonal.
%
% coeff is relative regularization w.r.t. the norm, e.g., 1e-6.

normX = normest(X, 1e-4);
n = length(X);
I = diag(sparse(ones(n,1)));

success = false;
times = 0;
%coeff = 1e-10;
reg = coeff*normX;
while ~success
  try
    Xreg = X + reg*I;
    U = chol(Xreg, varargin{:});
    success = true;
  catch
    reg = 100 * reg;
% $$$     if times == 0
% $$$       reg = 1e-16*normX;
% $$$     end
    times = times + 1;
% $$$     warning(sprintf(['Cholesky decomposition failed %d-th time, increasing' ...
% $$$                      ' regularization'], times))
    warning(sprintf(['Cholesky decomposition failed %d-th time, increasing ' ...
                     'regularization, coeff=%.1e, norm=%.1e'], times, ...
                    coeff, normX))
 end
end

%U = cholmod(sparse(X));
%U = U*S'; % re-permutate

% $$$ [L,D] = ldl(X);
% $$$ D = diag(D);
% $$$ D(D<0) = 0;
% $$$ U = diag(sqrt(D)) * L';

% $$$ choltol = 1e-6; % TODO: How to choose this??
% $$$ U = chol(X + diag(choltol*ones(length(X),1)));

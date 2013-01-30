% SPARSITY - Evaluates the sparsity of a matrix.
%
%   D = SPARSITY(X)
%
% D is the relative number of non-zero elements in the matrix X. For
% instance, D=0 for a matrix of zeros and D=1 for a matrix of no zeros.

function d = sparsity(X)

if issparse(X)
  d = nnz(X) / numel(X);
else
  d = sum(~isnan(X(:))) / numel(X);
end
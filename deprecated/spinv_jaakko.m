function Z = spinv(X)
% Calculates the sparse inverse of X, that is, only those elements (i,j)
% of the inverse that are non-zero in X. Note, that the real inverse can
% be full! ALSO, X should be symmetric!! What about positive definite?

warning('Use Vanhatalo''s sinv.m');

if ~issparse(X)
  error('X should be sparse');
end

[n,n] = size(X);
%Z = spones(X);%spalloc(M,N, nnz(X));
Z = X; % copy sparsity

L = chol(X, 'lower');
D = diag(L);
%I = spones(L);
%N =

[J,I] = find(tril(X));

% $$$ for i=N:-1:1
% $$$   for j=N:-1:i
for k=length(I):-1:1
  i = I(k);
  j = J(k);
    Z(i,j) = (i==j)/D(i)^2 - 1/D(i) * L((i+1):n,i)'*Z((i+1):n,j);
    Z(j,i) = Z(i,j);
end
% $$$   end
% $$$ end

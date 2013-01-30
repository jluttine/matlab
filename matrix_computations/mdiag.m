% D = MDIAG(X)
%
% Returns the diagonals of the matrices in X.
%
% X is M x N x K matrix.
%
% D is min(M,N) x K matrix having the diagonals of X(:,:,k) k=1,...,K as
% column vectors.

% Copyright (c) 2010 Jaakko Luttinen

function D = mdiag(X)

[M,N,K] = size(X);

D = zeros(min(M,N),K);
for k=1:K
  D(:,k) = diag(X(:,:,k));
end

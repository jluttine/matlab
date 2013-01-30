function X = mtraceprod(A,B)
% X = MTRACEPROD(A, B)
%
% Evaluates trace(A*B') efficiently.
%
% X = MTRACEPROD(A)
% 
% Evaluates trace(A*A') efficiently.
%
% If A and/or B has third dimension, the third dimension of A results in
% rows to X and the third dimension of B columns to X:
%
% X(I,J) = trace(A(:,:,I)*B(:,:,J)')

%warning('Deprecated.');

if nargin < 2
  [M,N,D] = size(A);
  X = reshape(A,[M*N,D])' * reshape(A,[M*N,D]);
else
  [M,N,Da] = size(A);
  [M,N,Db] = size(B);
  X = reshape(A,[M*N,Da])' * reshape(B,[M*N,Db]);
end
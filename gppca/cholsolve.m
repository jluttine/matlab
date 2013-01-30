
function X = cholsolve(A,B)
% Solves (A'*A)*X=B with respect to X. 
% A should be an upper triangular matrix (Cholesky decomposition).

opts.UT = true;
opts.TRANSA = true;
Z = linsolve(A, B, opts);
opts.TRANSA = false;
X = linsolve(A, Z, opts);


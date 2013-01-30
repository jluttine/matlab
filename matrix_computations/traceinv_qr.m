function y = traceinv_qr(Q,R)

A = linsolve_triu(R,speye(size(R)));
y = traceprod(A,Q);
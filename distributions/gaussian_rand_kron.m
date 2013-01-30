
function y = gaussian_rand_kron(L1,L2)
n1 = size(L1,2);
n2 = size(L2,2);
y = kronprod(L1,L2,randn(n2,n1));
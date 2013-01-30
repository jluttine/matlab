function test_conjgradmv()

N = 100;
x = rand(N,1);
A = sparse(double(rand(N) < 0.1));
K = A*A' + speye(N);

ind = 1:2:N;

xmv = x(ind);
Kmv = K(ind,ind);
ymv = Kmv*xmv;
y = zeros(N,1);
y(ind) = ymv;
I = false(N,1);
I(ind) = true;

x_pcg = pcg(Kmv,ymv,1e-6,100);
x_cg = conjgradmv(K,y,I,'tol',1e-6,'maxiter',100,'verbose',true);

norm(x_pcg-x_cg(ind))

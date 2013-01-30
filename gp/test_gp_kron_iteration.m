
function test_gp_kron_iteration

%
% DATA
%

logtheta1 = log(5);
N1 = 40;
x1 = 1:N1;
K1 = gp_cov_pp(logtheta1, x1, x1);

logtheta2 = log(5);
N2 = 45;
x2 = 1:N2;
K2 = gp_cov_pp(logtheta2, x2, x2);

s = 1e-0;
[X2,X1] = meshgrid(x2,x1);
freq1 = 0.03*2*pi;
freq2 = 0.07*2*pi;
Y_noiseless = lchol(K2) * randn(N2,N1) * lchol(K1)';
Y = Y_noiseless + s*randn(N2,N1);

%
% Comparison
%

K = kron(K1,K2) + s^2*speye(N1*N2);
tic
LD = ldlchol(K);
X_real = reshape(linsolve_ldlchol(LD, Y(:)), [N2,N1]);
toc

e = norm(Y);
V = Y/e;
for i=2:10
  W = K2*V*K1 + V;
  for j=1:k-1
    H(j,k-1) = traceprod(V,W);
    W = W - V*H(j,k-1);
  end
  H(k,k-1) = norm(W);
  V = W / H(k,k-1);
end


% $$$ tic
% $$$ LD1 = ldlchol(K1);
% $$$ LD2 = ldlchol(K2);
% $$$ X = zeros(N2,N1);
% $$$ for i=1:10
% $$$   Y1 = Y - s^2*X;
% $$$   Y2 = Y - kronprod(K1,K2,X);
% $$$   
% $$$   X1 = linsolve_kron_ldlchol(LD1,LD2,Y1);
% $$$   X2 = s^(-2)*Y2;
% $$$   
% $$$   X = 0.5*(X1+X2);
% $$$ end
% $$$ toc

norm(X-X_real)
function test_comparison_to_scipy

rand('state', 10)
randn('state', 10)

%M = 100
%N = 200
%D = 5
M = 200
N = 2000
D = 20

W = randn(D,M);
X = randn(D,N);
F = W'*X;
Y = F + 1*randn(M,N);

save('/home/jluttine/matlab/fa/data_pca_01_f.txt', 'F', '-ascii')
save('/home/jluttine/matlab/fa/data_pca_01_y.txt', 'Y', '-ascii')
save('/home/jluttine/matlab/fa/data_pca_01_d.txt', 'D', '-ascii')

Y(:,21:40) = NaN;

W_module = factor_module_iid(D,M);
X_module = factor_module_iid(D,N);
noise_module = noise_module_isotropic(M,N, ...
                                      'prior', struct('a_tau', 1e-5, ...
                                                  'b_tau', 1e-5));

Q = vbfa(D, Y, W_module, X_module, noise_module, ...
         'rotate', false, ...
         'maxiter', 200);

wx = Q.W'*Q.X;
WW = bsxfun(@times, reshape(Q.W,[D,1,M]), reshape(Q.W,[1,D,M])) + Q.CovW;
XX = bsxfun(@times, reshape(Q.X,[D,1,N]), reshape(Q.X,[1,D,N])) + Q.CovX;
wwxx = reshape(WW, [D*D,M])' * reshape(XX, [D*D,N]);
err_wx = 2*sqrt(wwxx - wx.^2);

tserrorplot(wx', err_wx')
addtsplot(W'*X, 'g')

%Q.Tau

%Q.CovX
%Q.CovW
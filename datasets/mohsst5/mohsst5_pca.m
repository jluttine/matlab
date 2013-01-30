function Q = mohsst5_pca(data, D, maxiter)

if nargin < 3
  maxiter = 10;
end

Q = pca_full(data.observations,D,'maxiters',maxiter,'rotate2pca',true, ...
                 'algorithm','ppca');

Yh_ppca = Q.A * Q.S + repmat(Q.Mu,1,N);

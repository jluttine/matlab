
function test_preconditioner(type)

N = 1500;
x = 1:N;

logtheta = log(20);
K = gp_cov_pp(logtheta, x, x);

imv = (zeros(N,1) ~= 0);
%imv(2:3:N) = true;

LD = ldlchol(K);
%[L,D] = ldlsplit(LD);
L = ldlchol2chol(LD);

I = speye(N);

s = 1e-4;
y = L*randn(N,1);
ynoise = y + s*randn(N,1);


%inv_Dmv = inv(D(imv,imv));
%x = I(~imv,:)*ldlsolve(LD,I(:,~imv)*yobs) ...
%    - I(~imv,:)*ldlsolve(LD, I(:,imv) * (inv_Dmv \ (I(imv,:)*ldlsolve(LD, ...
%                                                  I(:,~imv)*yobs))))

% Different preconditioners:

Cov = K + s^2*I;
LD_Cov = ldlchol(Cov);

f1 = @(x) ldlsolve(LD_Cov, x);

get_f = @(invK) ( @(x) f1(x) ...
                  - ldlsolve(LD_Cov, I(:,imv)*invK((I(imv,:)*ldlsolve(LD_Cov,x)))) );

maxiter = 1000;
switch type
 case 0
  f = [];
 case 1
  f = f1;
 case 2
  g = @(x) Cov(imv,imv)*x;
  f = get_f(g);
 case 3
  g = @(x) diag(diag(Cov(imv,imv)))*x;
  f = get_f(g);
 case 4
  invCov = inv(full(Cov));
  g = @(x) diag(diag( invCov(imv,imv) )) \ x;
  f = get_f(g);
 case 5
  % Exact:
  %g = @(x) Cov(imv,imv)*x - Cov(imv,~imv)*(Cov(~imv,~imv) \ (Cov(~imv,imv)*x));
  % Something between?:
  %g = @(x) Cov(imv,imv)*x - Cov(imv,~imv)*I(~imv,:)*(ldlsolve(LD_Cov,(I(:,~imv)*Cov(~imv,imv)*x)));
  % Diagonalized:
  %g = @(x) Cov(imv,imv)*x - Cov(imv,~imv)*(diag(diag(Cov(~imv,~imv))) \ (Cov(~imv,imv)*x));
  f = get_f(g);
 case 6
  % Neumann
  h = @(x) I(~imv,:)*(Cov*(I(:,~imv)*x));
  f = @(x) I(:,~imv)*linsolve_neumann(h, I(~imv,:)*x, normest(Cov), 100);
  maxiter = 10;
 case 7
  % reduced rank with svd
  [U,S] = svds(Cov(~imv,~imv),30);
  size(Cov(~imv,~imv))
  d = diag(Cov(~imv,~imv)) - dot(U,U,2);
  mean(d)
  f = @(x) I(:,~imv) * ((U*U'+spdiag(d))\(I(~imv,:)*x));
 case 8
  % reduced rank with inducing inputs
  u = x(:,1:5:end);
  K_u = gp_cov_pp(logtheta, u, u);
  K_xu = gp_cov_pp(logtheta, x, u);
  D = spdiag( diag(K) - diag(K_xu*(K_u\K_xu')) ) + s^2*speye(N);
  D = D(~imv,~imv);
  LD_u = ldlchol(K_u+K_xu'*I(:,~imv)*inv(D)*I(~imv,:)*K_xu);
  f = @(x) I(:,~imv)*(D\(I(~imv,:)*x) - D\(I(~imv,:)*K_xu* ...
                                           linsolve_ldlchol(LD_u, ...
                                                    K_xu'*I(:,~imv)* ...
                                                    inv(D)*I(~imv,:)* ...
                                                    x)));
 case 9
  % block diagonal approximation
  % O(M^2*N)
  blksize = 20;
  Cov_blk = spalloc(N,N,blksize*N);
  for i=1:blksize:N
    j = i+blksize-1;
    if j <= N
      Cov_blk(i:j,i:j) = Cov(i:j,i:j);
    else
      Cov_blk(i:end,i:end) = Cov(i:end,i:end);
    end
  end
  LD_blk = ldlchol(I(~imv,:)*(Cov_blk+s^2*speye(N))*I(:,~imv));
  f = @(x) I(:,~imv)*linsolve_ldlchol(LD_blk,I(~imv,:)*x);
  
 case 10
  % reduced cholesky
  L = ldlchol2chol(LD);
  Lo = L(~imv,~imv);
  Lm = L(~imv,imv);
  LDo = LD(~imv,~imv);
  LDm = L(~imv,imv);
  %clf, imagesc(Lo), pause, imagesc(inv(Lo)); pause
  
  %Ko = Lo*Lo';
  %clf, imagesc(Ko), return
  
  %clf, imagesc(linsolve_chol(Lo,Lm,'lower')); return
  Z = eye(sum(imv)) ;% + Lm'*(linsolve_chol(Lo,Lm,'lower'));
  %Z = eye(sum(imv)) + Lm'*(linsolve_ldlchol(LDo,Lm));
  %clf, imagesc(Z), return
  L_z = chol(Z, 'lower');
  f = @(x) I(:,~imv)*(linsolve_chol(Lo,I(~imv,:)*x, 'lower') - ...
                      linsolve_chol(Lo,Lm*linsolve(L_z,Lm* ...
                                                   linsolve_chol(Lo,I(~imv,:)*x, ...
                                                    'lower'), ...
                                                   'lower'), ...
                                    'lower') );

  % ignore removed stuff
  f = @(x) I(:,~imv)*linsolve_chol(Lo,I(~imv,:)*x, 'lower');
  f = @(x) I(:,~imv)*linsolve_ldlchol(LDo,I(~imv,:)*x);

 case 11
  % evaluate diagonal (or also some other) elements of the inverse matrix
  % brute force for testing..
  D = diag(diag( inv(Cov(~imv,~imv)) ));
  f = @(x) I(:,~imv)*D*I(~imv,:)*x;
  
  invCov = inv(Cov(~imv,~imv));
  %clf, subplot(2,1,1), imagesc(invCov), subplot(2,1,2), imagesc(Cov),return
  nobs = sum(~imv);
  l = 1;
  M = spdiags(ones(nobs,2*l+1), -l:l, nobs, nobs) ~= 0;
  invCov = M.*invCov; % mask
  f = @(x) I(:,~imv)*invCov*I(~imv,:)*x;
  
  
 case 12
  % kiintopisteiteraatio..
  
  
  A = Cov(~imv,~imv);
  B = Cov(~imv,imv);
  C = Cov(imv,~imv);
  D = Cov(imv,imv);
  invCov = inv(Cov);
  Ai = invCov(~imv,~imv);
  Bi = invCov(~imv,imv);
  Ci = invCov(imv,~imv);
  Di = invCov(imv,imv);
  
  z = randn(sum(~imv),1);
  z = ones(sum(~imv),1);
  v = zeros(sum(~imv),1);%z
  r = 1e-0;
  v = r*(A\z + 0.01*randn(size(z)));
  for i=1:10
    v = r*( Ai*z - Bi*(D*(Ci*z)) - Bi*(C*(Ai*z)) + Bi*(C*v) );
  end
  v/r
  A\z
  return
  
 case 13
  
  % combined block-diagonal + inducing inputs
  
  % Ni = inducing inputs in a block
  % BLKSIZE = size of the blocks
  % you want at least that Ni^2 <= BLKSIZE (doesn't use more memory than
  % data)
  
  % Problem: the number of inducing inputs needs to be very large.. :(
  
  % block diagonal approximation
  blksize = N;
  Nblk = N/blksize;
  Ni = 50;%floor(sqrt(blksize));
  L_ublks = zeros([Ni,Ni,Nblk]);
  %Cov_blk = spalloc(N,N,blksize*N);
  ui = [];
  Kblk = spalloc(N,N,N*blksize);
  Iblk = speye(blksize);
  for n=1:Nblk
    i0 = (n-1)*blksize + 1;
    i1 = min(n*blksize,N);
    i = round(linspace(i0-1,i1+1,Ni+2));
    i = i(2:(end-1));
    K_u = full(K(i,i));
    K_ux = K(i,i0:i1);
    KUX{n} = K_ux;
    K_approx = K_ux'*linsolve_cov(K_u,K_ux);
    %full(K(i0:i1,i0:i1))
    D = spdiag(diag(K(i0:i1,i0:i1)) - diag(K_approx)) + s^2*Iblk;
    D = D(~imv(i0:i1),~imv(i0:i1));
    invD{n} = inv(D);
    Cov_u = K_u + K_ux*Iblk(:,~imv(i0:i1))*invD{n}*Iblk(~imv(i0:i1),:)*K_ux';
    %subplot(2,1,1), imagesc(K_approx), subplot(2,1,2), imagesc(K(i,i)), ...
    %        return
    L_covu{n} = sparse(chol(Cov_u, 'lower'));
    ui = [ui,i];
    
    % true block diag for debugging
    Kblk(i0:i1,i0:i1) = K(i0:i1,i0:i1);
  end
  L_covu = blkdiag(L_covu{:});
  invD = blkdiag(invD{:});
  K_ux = blkdiag(KUX{:});
  
  % These should be identical:
% $$$   clf
% $$$   subplot(3,1,1), imagesc(invD-invD*I(~imv,:)*K_ux'*inv(L_covu*L_covu')*K_ux*I(:,~imv)*invD)
% $$$   subplot(3,1,2), imagesc(inv(Kblk(~imv,~imv)+s^2*speye(sum(~imv))))
% $$$   subplot(3,1,3), imagesc(inv(K(~imv,~imv)+s^2*speye(sum(~imv))))
% $$$   return
  
%  LD_blk = ldlchol(I(~imv,:)*(Cov_blk+s^2*speye(N))*I(:,~imv));
  f = @(x) I(:,~imv)*(invD*(I(~imv,:)*x) - invD*(I(~imv,:)*(K_ux'* ...
                                                    (linsolve_chol(L_covu',K_ux*(I(:,~imv)*(invD*I(~imv,:)*x)))))));
  
 case 14
  % This gives exact..
  invCov = inv(full(Cov));
  g = @(x) invCov(imv,imv) \ x;
  f = get_f(g);
end

%K_joint = [K, K; K, K+s^2*I];
x = gaussian_rand_pcg(ynoise, Cov, K, L*randn(N,1), s*randn(N,1), f, ...
                      I(~imv,:), 'maxiter', maxiter);

yobs = ynoise;
yobs(imv) = nan;

clf
plot(y, 'k')
hold on
plot(yobs, 'rx')
plot(x, 'g');


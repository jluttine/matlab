function test_lchol_grad()

N = 100;
%x = randn(N,1)';%1:N;
x = 1:N;
D = sqrt(sq_dist(x));
covfunc = gp_cov_scale(gp_cov_pp(D,1));
[K,dK] = covfunc([1 30]);
%[K,dK] = covfunc(0.02);
% $$$ figure
% $$$ imagesc(K)
% $$$ figure
% $$$ imagesc(D)
% $$$ return

% return

LD = ldlchol(K);
% $$$ [LD,p,q] = ldlchol(K);
I = find(LD);
%J = find(tril(dK{1}));
% Z = sparse(I,J,1);% - spones(LD);
L = ldlchol2lchol(LD);

% $$$ %[S Z] = sparse2(ldlchol(K));
% $$$ figure
% $$$ spy(Z)
% $$$ figure
% $$$ spy(L)
% $$$ % $$$ figure
% $$$ % $$$ imagesc(LD)
% $$$ return

% $$$ [~,~,~,~,Z] = symbfact(K);
% $$$ figure
% $$$ imagesc(Z');

% $$$ return
% $$$ 
t = cputime();
L = lchol(K);
fprintf('Time for Cholesky: %.3f s\n', cputime()-t);

% $$$ dL = cell(length(dK),1);
% $$$ for m=1:length(dK)
% $$$   I = tril(ones(size(L))) - 0.5*eye(size(L));
% $$$   dL{m} = 0.5 * (L \ (I.*dK{m}));
% $$$ end

mex lchol_grad_sparse.c
t = cputime();
[L, dL] = lchol_grad(K, dK);
%dL = lchol_grad_sparse(dK{1}, L, I);
fprintf('Time for Cholesky gradient: %.3f s\n', cputime()-t);

for m=1:length(dK)
  err = dK{m} - dL{m}*L' - L*dL{m}';
  normest(err)
  full(max(abs(err(:))))
end

return

% $$$ figure
% $$$ spy(LD)
% $$$ 
% $$$ figure
% $$$ spy(dK{1})
% $$$ 
for m=1:length(dL)
  figure
  imagesc(dL{m})
  colormap_redblue();
  colormap_centered();

  rec = dL{m}*L' + L*dL{m}';
  figure
  imagesc(rec)
  colormap_redblue();
  colormap_centered();
end

% $$$ truth = dK{1};
% $$$ figure
% $$$ imagesc(truth)
% $$$ colormap_redblue();
% $$$ colormap_centered();

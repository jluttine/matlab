% [L, dL] = lchol_grad(K, dK)
%
% Returns the Cholesky decomposition of a pos.def.symm. matrix and the
% derivatives of the factor. 

% TODO: This is not optimized! Just a quick and dirty solution. Also, the
% mex file is not optimized..
function [L, dL] = lchol_grad(K, dK)

N = size(K,1);

% Find symbolic structure
[cols,~,~,~,Z] = symbfact(K,'sym','lower');
cols = [0, cumsum(cols)];

I = find(Z);

% Transpose
% I = sort(mod(I-1,N)*N + ceil(I/N));

% Form a matrix of all data
if iscell(dK)
  M = numel(dK);
else
  M = 1;
end
V = zeros(length(I), 1+1+M);
V(:,1) = I;
V(:,2) = K(I);
if iscell(dK)
  for m=1:M
    V(:,2+m) = dK{m}(I);
  end
else
  V(:,3) = dK(I);
end

% $$$ figure
% $$$ imagesc(V(:,2:end));

% Solve using a MEX function
try
  Y = lchol_grad_sparse(V, cols, N);
catch
  L = NaN;
  dL = NaN;
  return;
end

i = mod(I-1,N)+1;
j = floor((I-1)/N)+1;
L = sparse(i,j,Y(:,2));
if iscell(dK)
  dL = cell(M,1);
  for m=1:M
    dL{m} = sparse(i,j,Y(:,2+m));
  end
else
  dL = sparse(i,j,Y(:,3));
end

% $$$ figure
% $$$ imagesc(Y(:,2:end));

% $$$ I = find(Z');
% $$$ figure
% $$$ plot(I-It)

%Z = Z';
%I = find(Z);



return

N = size(L,1);

dV = full(dK(I)); % ???

%tmp = [];

ind = 1;
for k=1:N

  i = mod(I(ind)-1,N)+1;
  j = ceil(I(ind)/N);
  if i~=j || i~=k
    error('WHAAAAT');
  end

  % Pivot
  dV(ind) =  0.5 * (dV(ind)/L(k,k));
  ind_pivot = ind;
  
  % Lead column
  ind = ind_pivot + 1;
  while ind <= length(I) && I(ind) <= k*N
    j = mod(ind-1,N) + 1;
    dV(ind) = (dV(ind) - L(j,k)*dV(ind_pivot)) / L(k,k);
    ind = ind + 1;
  end

  % Rows
  ind = ind_pivot + 1;
  while ind <= length(I) && I(ind) <= k*N
    i = mod(I(ind)-1,N) + 1; % column
    jnd = ind;
    while jnd <= length(I) && I(jnd) <= k*N
      j = mod(I(jnd)-1,N) + 1; % row
      %j = ceil(I(jnd)/N); % column
      l = find(I==(i-1)*N+j, 1);
      if isempty(l)
        I(:)
        (i-1)*N+j
        i
        j
        error('WHAAAT');
      end
      if l <= jnd
        error('whatataaaa');
      end
      dV(l) = dV(l) - dV(ind)*L(j,k) - dV(jnd)*L(i,k);
      jnd = jnd + 1;
    end
    ind = ind + 1;
  end
  
%  tmp = [tmp, dV];
  
end

%tmp

i = mod(I-1,N)+1;
j = ceil(I/N);
dL = sparse(i,j,dV,N,N);

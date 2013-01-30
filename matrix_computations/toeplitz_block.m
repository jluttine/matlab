% TOEPLITZ_BLOCK - Constructs a block-Toeplitz matrix.
%
%   Y = TOEPLITZ_BLOCK(X)
%
% X = [X_1, X_2, X_3, ..., X_N]
%
% Y = [X_1 X_2 X_3 ... X_N
%      X_2 X_1 X_2 ... ...
%      X_3 X_2 X_1 ... ...
%      ... ... ... ... ...
%      X_N ... X_3 X_2 X_1]


function Y = toeplitz_block(X)

[M,N] = size(X);

if mod(N,M) ~= 0
  error('The number of columns must be a multiple of the number of rows.')
end

D = N/M;

if issparse(X)
  t = cputime();
  
  nonzeros = zeros(D,1);
  i = cell(D,1);
  j = cell(D,1);
  v = cell(D,1);
  [i{1},j{1},v{1}] = find(X(:,1:M));
  nonzeros(1) = D*length(v{1});
  for d=2:D
    k = (d-1)*M+1;
    l = d*M;
    [i{d},j{d},v{d}] = find(X(:,k:l));
    nonzeros(d) = 2 * (D-d+1) * length(v{d});
  end
  nzs = sum(nonzeros);
  
  
  %time = cputime() - t, t = cputime();
  
  I = zeros(nzs,1);
  J = zeros(nzs,1);
  values = zeros(nzs,1);
  
  %time = cputime() - t, t = cputime();

  z = 1;
  for d2=1:D
    for d1=1:D
      ind = abs(d1-d2) + 1;
      jnd = z:(z+length(v{ind})-1);
      I(jnd) = i{ind} + (d1-1)*M;
      J(jnd) = j{ind} + (d2-1)*M;
      values(jnd) = v{ind};
      z = z + length(v{ind});
    end
  end

  %time = cputime() - t, t = cputime();
  
  Y = sparse(I,J,values,N,N,nzs);
  
  %time = cputime() - t, t = cputime();
  

else
  
  Y = zeros(N,N);
  for i=1:D
    for j=1:D
      ind = abs(i-j) + 1;
      k = (i-1)*M+1;
      l = (j-1)*M+1;
      m = (ind-1)*M+1;
      Y(k:(k+M-1),l:(l+M-1)) = X(:,m:(m+M-1));
    end
  end
  
end


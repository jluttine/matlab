
function test_sparse_chol
warning('This function is deprecated')


N = 2000;
x = 1:N;
K = gp_cov_pp(log(100), x, x);

imagesc(K)

disp('Cholesky')
tic
  L = chol(K,'lower');
  U = L';
toc

tic
  [LD,p,q] = ldlchol(K);
toc
nnzld = nnz(LD)
tic
  LD = ldlchol(K);
toc
nnzld = nnz(LD)

y = randn(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 1000;

disp('Solving vector simply and using Cholesky')

% $$$ t = cputime;
% $$$ for m=1:M
% $$$   z = K \ y;
% $$$ end
% $$$ fprintf('Elapsed time is %.8f\n', (cputime-t)/M);
  
t = cputime;
for m=1:M
  z = L\y;
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  z = L'\z;
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  z = L\y;
  z = U\z;
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  z = L\y;
  z = L'\z;
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  z = ldlsolve(LD,y);
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 5;

Y = randn(N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Solving matrix simply and using Cholesky')

t = cputime;
for m=1:M
  Z = K \ Y;
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  % For some reason it is faster when dividing into two separate rows..?
  Z = U\(L\Y);
  %Z = L\Y;
  %Z = U\Z;
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  Z = ldlsolve(LD,Y); % is the fastest!
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

% $$$ invK = inv(K);
% $$$ t = cputime;
% $$$ for m=1:M
% $$$   Z = invK*Y;
% $$$ end
% $$$ fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Multiplying from left.')

t = cputime;
for m=1:M
  K*Y; % a bit slower
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  K'*Y; % a bit faster
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Multiplying from right.')

M = 10;

t = cputime;
for m=1:M
  Y*K'; % a bit slower
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  Y*K; % a bit faster
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Multiplying a vector from left.')

M = 1000;

t = cputime;
for m=1:M
  K*y; % slower
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);

t = cputime;
for m=1:M
  K'*y; % faster
end
fprintf('Elapsed time is %.8f\n', (cputime-t)/M);


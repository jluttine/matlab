function Q = alex2jaakko(A)
Q.W = A.A;
Q.X = A.S;
Q.mu = A.Mu;

if ~isempty(A.Av)
  D = rows(A.Av{1});
  M = length(A.Av);
  Q.CovW = zeros([D,D,M]);
  for m=1:M
    Q.CovW(:,:,m) = A.Av{m};
  end
end
if ~isempty(A.Sv)
  D = rows(A.Sv{1});
  N = length(A.Sv);
  Q.CovX = zeros([D,D,N]);
  for n=1:N
    Q.CovW(:,:,n) = A.Sv{n};
  end
end

if isempty(A.Muv)
  Q.v_mu = 0;
else
  Q.v_mu = A.Muv;
end

function [W,X,S,mu,v,tau] = rppca(Y, d, varargin)
% Implementation of robust probabilistic PCA.
% See Archambeau, Delannay, Verleysen: "Robust Probabilistic
% Projections", 2006.

opts = struct( ...
    'init',          [],...
    'prior', [], ...
    'hiermean', false, ...
    'rotate', true, ...
    'commonnoise', true, ...
    'startrotate', 1, ...
    'startupdatehyper', 1, ...
    'autosavetime', 0,...
    'autosavefile', 'vbpcamv_autosave',...
    'testset', [], ...
    'fixw', false, ...
    'maxiters', 1000, ...
    'convergence', eps);


[D,n] = size(Y);

% initializing guesses
mu = zeros(D,1);%mean(Y,2);
W = orth(randn(D,d));%rand(D,d);
tau = 1;
v = 10;

% initialize sizes
u = zeros(1,n);
logu = zeros(1,n);
X = zeros(d,n);
S = zeros(d,d,n);
ID = eye(D);
Id = eye(d);

for i=1:100
  
  % E step
  A = inv(W*W' + ID/tau);
  B = tau * W'*W + Id;
  for j=1:n
    u(j) = (D+v) / ((Y(:,j)-mu)'*A*(Y(:,j)-mu)+v);
    logu(j) = psi((D+v)/2) - log(((Y(:,j)-mu)'*A*(Y(:,j)-mu)+v)/2);
    X(:,j) = tau * inv(B) * W' * (Y(:,j) - mu);
    S(:,:,j) = inv(B) + u(j)*X(:,j)*X(:,j)';
%    Sv{j} = inv(B) / u(j);
  end
  
%  [dmu,W,X,Sv] = RotateToPCA(W,X,tau*inv(B),[]);
% $$$   [dmu,W,X,Sv] = RotateToPCA(W,X,Sv,[]);
% $$$   for j=1:n
% $$$     S(:,:,j) = u(j)*Sv{j} + u(j)*X(:,j)*X(:,j)';
% $$$   end
  
  % M step
  mu = ((Y-W*X)*u') / sum(u);
  Mu = repmat(mu,1,n);
  W = zeros(size(W));
  for j=1:n
    W = W + u(j)*(Y(:,j)-mu)*X(:,j)';
  end
  W = W * inv(sum(S,3));
%  W = ((U.*(Y-Mu))*X') * inv(sum(S,3));
  SWW = 0;
  for j=1:n
    SWW = SWW + trace(S(:,:,j) * W'*W);
  end
  c = 0;
  for j = 1:n
    c = c - 2*u(j)*(Y(:,j)-mu)'*W*X(:,j);
  end
  N = n * D;
  tau = 1/(sum(((Y-Mu).^2)*u')/N + c/N + SWW/N);
  c = mean(logu-u);
  
  v = exp(fminsearch(@(v) ((1+log(exp(v)/2)-psi(exp(v)/2)+c)^2), ...
                     log(v)));
  
end


function [ dMu, A, S, Sv ] = RotateToPCA( A, S, Sv, Isv, obscombj );

n2 = size(S,2);

mS = mean(S,2);
dMu = A*mS;
S = S - repmat(mS,1,n2);

covS = S*S';
if isempty(Isv)
  %covS = n2*Sv;
    for j = 1:n2
        covS = covS + Sv{j};
    end
else
    nobscomb = length(obscombj);
    for j = 1:nobscomb
        covS = covS + ( length(obscombj{j})*Sv{j} );
    end
end

covS = covS / n2;
[VS,D] = eig(covS);
A = A*VS*sqrt(D);

[ A, DA, VA ] = svd(A);
A = A*DA;
R = VA'*diag(1./sqrt(diag(D)))*VS';

S = R*S;
%Sv = R*Sv*R';
for j = 1:length(Sv)
    Sv{j} = R*Sv{j}*R';
end

function [W,CovW,X,CovX,R] = rotate_to_pca(W,CovW,X,CovX,weightsWW)
% [W,CovW,X,CovX,R] = rotate_to_pca(W,CovW,X,CovX,weightsW)

% Dimensionalities
[M,D] = size(W);
[D,N] = size(X);

% Convert cell covariance matrices to arrays
cellW = false;
if iscell(CovW)
  CovW = covcell_to_covarray(CovW);
  cellW = true;
end

cellX = false;
if iscell(CovX)
  CovX = covcell_to_covarray(CovX);
  cellX = true;
end

if nargin < 5 || isempty(weightsWW)
  weightsWW = ones(M,1);
else
  weightsWW = weightsWW(:) .* ones(M,1);
end

%% Find mixing rotation Xpca = R * Xgp

R = 1;
    
% 1) Whiten <XX'>

XX = X * X' + sum(CovX,3);
[V,A] = svd(XX/N);
R = diag(sqrt(1./diag(A))) * V';
% Rotate W
W = W / R;
for i=1:rows(W)
  CovW(:,:,i) = R' \ CovW(:,:,i) / R;
end

% 2) Orthogonalise weighted <W'W>

% Evaluate weighted second moment
WW = W' * diag(weightsWW) * W;
for m=1:M
  WW = WW + weightsWW(m) * CovW(:,:,m);
end

% Diagonalise it
[V,D] = svd(WW);
R = V' * R;

% Rotate W
W = W * V;
for m=1:M
  CovW(:,:,m) = V' * CovW(:,:,m) * V;
end

% $$$ % Debug the rotation!
% $$$ weights2 = cosd(data.coordinates(2,:));
% $$$ WW = W' * diag(weights2) * W;
% $$$ for m=1:M
% $$$   WW = WW + weights2(m) * CovW(:,:,m);
% $$$ end
% $$$ XX = R * (X*X'+sum(CovX,3)) * R';
% $$$ WW(1:10,1:10)
% $$$ XX(1:10,1:10) / size(X,2)


% Rotate X to PCA
X = R * X;
for n=1:N
  CovX(:,:,n) = R * CovX(:,:,n) * R';
end

% Convert back to cells
if cellW
  CovW = covarray_to_covcell(CovW);
end
if cellX
  CovX = covarray_to_covcell(CovX);
end
function [W,CovW,X,CovX] = mohsst5_rotate_to_pca(W,CovW,X,CovX)

% Remove temporal mean
mu_X = mean(X,2);
X = bsxfun(@minus, X, mu_X);

% Compute <XX>
XX = X*X';
if ndims(CovX) == 2 && all(size(CovX)==size(X))
  XX = XX + diag(sum(CovX,2));
else
  XX = XX + sum(CovX,3);
end

% Compute weighted <WW>
w = mohsst5_weights();
w = mohsst5_remove_land(w);
WW = W*diag(w)*W';
if ndims(CovW) == 2 && all(size(CovW)==size(W))
  WW = WW + diag(wsum(CovW,w,2));
else
  WW = WW + wsum(CovW,w,3);
end

% Find the SVD of the rotation:
% R = Ur * Dr * Vr'

% Whiten <XX>
[Ux,Dx] = svd(XX/size(X,2));
Vr = Ux;
Dr = diag(diag(Dx).^(-0.5));

% Orthogonalize <WW>
[Uw,~] = svd(inv(Dr)*Vr'*WW*Vr*inv(Dr));
Ur = Uw';

R = Ur*Dr*Vr';

% $$$ figure
% $$$ imagesc(XX)
% $$$ figure
% $$$ imagesc(R*XX*R')
% $$$ figure
% $$$ imagesc(WW)
% $$$ figure
% $$$ imagesc(R'\WW/R)
% $$$ return

X = R*X;
W = R'\W;

if ndims(CovX) == 2 && all(size(CovX)==size(X))
  CovX = R.^2*CovX;
else
  for n=1:size(CovX,3)
    CovX(:,:,n) = R*CovX(:,:,n)*R';
  end
end
  
if ndims(CovW) == 2 && all(size(CovW)==size(W))
  CovW = inv(R').^2*CovW;
else
  for n=1:size(CovW,3)
    CovW(:,:,n) = R'\CovW(:,:,n)/R;
  end
end

X = bsxfun(@plus, X, R*mu_X);

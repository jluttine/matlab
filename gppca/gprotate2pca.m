function [W,varW,X,varX] = gprotate2pca(W,varW,X,varX,weights)

[M,D] = size(W);
N = cols(X);

%WW = W'*W + diag(rowsum(varW));

% $$$ % Use ML zero mean (mu should be updated after this function!!) ..
% $$$ dmu = mean(X,2);

% $$$ % Move bias
% $$$ X = X - repmat(dmu,1,n);
% $$$ mu = mu + W*dmu;

Qx = 1;

%warning('DISCARDING VARIANCES!!')

% Whiten w.r.t. X
muX = mean(X,2);
X0 = bsxfun(@minus,X,muX);
XX = X0*X0' + diag(colsum(varX));
% $$$ if fixw
% $$$   [Vx,Dx,tmp] = svd(XX/(n-m)); % USE THIS IF FIXED w ??
% $$$ else
[Vx,Dx,tmp] = svd(XX/N);
% $$$ end
Qx = diag(1./sqrt(diag(Dx))) * Vx';
Qw = Vx*sqrt(Dx);
W = W * Qw;
for i=1:M
  varW(i,:) = diag(Qw'*diag(varW(i,:))*Qw);
end
X = Qx * X;
for j = 1:N
  varX(:,j) = diag(Qx*diag(varX(:,j))*Qx');
end

rotationX = Qx;

% $$$ % Check that XX is really whitened! (because of numerical issues!!)
% $$$ XX = X*X' + sum(CovX,3);
% $$$ if fixw
% $$$   [Vx,Dx,tmp] = svd(XX/(n-m)); % USE THIS IF FIXED w ??
% $$$ else
% $$$   [Vx,Dx,tmp] = svd(XX/n);
% $$$ end
% $$$ if Dx(1) > 1.1
% $$$   % Whiten w.r.t. X AGAIN
% $$$   Qx = diag(1./sqrt(diag(Dx))) * Vx';
% $$$   Qw = Vx*sqrt(Dx);
% $$$   W = W * Qw;
% $$$   for i=1:size(CovW,3)
% $$$     CovW(:,:,i) = Qw'*CovW(:,:,i)*Qw;
% $$$   end
% $$$   X = Qx * X;
% $$$   for j = 1:size(CovX,3)
% $$$     Sv{j} = Qx*Sv{j}*Qx';
% $$$     CovX(:,:,j) = Qx*CovX(:,:,j)*Qx';
% $$$   end
% $$$ end

Qx = 1;
% Diagonalize w.r.t. W
WW = W'*W + diag(rowsum(varW));
% Use weights!!
if nargin >= 5
  WW = W'*diag(weights(:))*W + diag(weights(:)' * varW);
%  WW = W'*diag(weights(:))*W + diag(rowsum(varW));
% $$$   w = sqrt(weights);
% $$$   WW = bsxfun(@times, WW, w(:));
% $$$   WW = bsxfun(@times, WW, w(:)');
end
%WW = W'*W + sum(CovW,3);
[Vw,Dw,tmp] = svd(WW);
%norms_of_W = sqrt(diag(Dw))
%[Dw,I] = sort(diag(Dw), 'descend');
%Vw = Vw(:,I);
Qx = Vw' * Qx;
Qw = Vw;
W = W * Qw;
for i=1:M
  varW(i,:) = diag(Qw'*diag(varW(i,:))*Qw);
end
X = Qx * X;
for j = 1:N
  varX(:,j) = diag(Qx*diag(varX(:,j))*Qx');
end

rotationX = Qx * rotationX;

% Rotate such that, the largest loadings are positive
for d=1:cols(W)
  if max(W(:,d)) < max(-W(:,d))
    W(:,d) = -W(:,d);
    X(d,:) = -X(d,:);
  end
end

% $$$ figure
% $$$ imagesc(abs(rotationX));
% $$$ 
% $$$ figure
% $$$ plot( sum(W.^2,1) )
% $$$ 
% $$$ figure
% $$$ plot( diag(Dw) )

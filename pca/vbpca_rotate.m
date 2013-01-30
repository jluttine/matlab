function [W,CovW,X,CovX] = vbpca_rotate(R, W, CovW, X, CovX)
% Rotate W*R and R\X

% Rotate W
W = W * R;
if ~isempty(CovW)
  if iscell(CovW)
    for i=1:length(CovW)
      CovW{i} = R'*CovW{i}*R;
    end
  else
    for i=1:size(CovW,3)
      CovW(:,:,i) = R'*CovW(:,:,i)*R;
    end
  end
end

% Rotate X
X = R \ X;
if ~isempty(CovX)
  if iscell(CovX)
    for j = 1:length(CovX)
      CovX{j} = (R \ (R\CovX{j})')';
    end
  else
    for j = 1:size(CovX,3)
      CovX(:,:,j) = (R \ (R\CovX(:,:,j))')';
    end
  end
end

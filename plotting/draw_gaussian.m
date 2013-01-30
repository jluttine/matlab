function draw_gaussian(mu,Cov,std,varargin)
if nargin < 3 || isempty(std)
  std = 1:3;
end

[U,S,V] = svd(Cov);
a = linspace(0,2*pi,100);%0:0.1:(2*pi);
X = [cos(a);sin(a)];

for ind=1:length(std)
  Y = std(ind) * U*sqrt(S)*X;
  Z1(:,ind) = mu(1) + Y(1,:);
  Z2(:,ind) = mu(2) + Y(2,:);
end

plot(Z1,Z2,varargin{:});

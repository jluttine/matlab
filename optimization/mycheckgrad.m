function d = mycheckgrad(f, x, e, varargin);
% d = mycheckgrad(f, x, e);

[y, dy] = feval(f, x, varargin{:});             % get the partial derivatives dy
n = length(x);

dh = zeros(n,1) ;
for j = 1:n
  dx = zeros(n,1);
  dx(j) = dx(j) + e;                               % perturb a single dimension
  y2 = feval(f, x+dx, varargin{:});
  y1 = feval(f, x-dx, varargin{:});
  dh(j) = (y2 - y1)/(2*e);
  fprintf('%+.4e \t %+.4e \t %+.4e \t (function, numerical, x)\n', dy(j), dh(j), x(j));
end

%disp([dy dh])                                          % print the two vectors
%[dy dh]
d = norm(dh-dy)/norm(dh+dy);       % return norm of diff divided by norm of sum

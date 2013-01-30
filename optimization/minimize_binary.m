function x = minimize_binary(f,df,x0)
% Binary search for minimizing a unimodal univariate function with a
% strictly monotonically increasing derivative.

step = 1;
dy0 = feval(df,x0);
s0 = sign(dy0);
s1 = s0;

% Find boundary [x0,x1] (or [x1,x0]) inside which the optimum lies.
while s1 ~= -s0
  x1 = x0;
  x0 = x0 - s0 * step;
  dy1 = dy0;
  dy0 = feval(df,x0);
  s0 = sign(dy0);
  step = 2*step;
end
[x0, x1]

x = fzero(

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Binary search within the boundary
for m=1:10
  x = 0.5*(x0+x1);
  dy = feval(df,x);
  if sign(dy) == s0
    x0 = x;
  else
    x1 = x;
  end
end
y1 = feval(f,x1);
y0 = feval(f,x0);

[x0, x1]

% Finally, try linear interpolation between the boundary points
x = (abs(dy1)*x0+abs(dy0)*x1)/(abs(dy1)+abs(dy0));
y = feval(f,x);

[x0, x1, x]

% Select the optimum
if y0 < y
  x = x0;
elseif y1 < y
  x = x1;
end


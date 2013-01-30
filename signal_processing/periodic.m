function y = periodic(x, x0, x1)

l = x1-x0;

y = rem(x-x0, l);
y(y<0) = y(y<0) + l;
y = y + x0;
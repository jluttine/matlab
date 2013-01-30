function [f,df] = tmpsin(x)
f = sum(sin(x));
df = cos(x);
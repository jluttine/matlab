function test_t_tml

n = 100;
t2 = trnd(10,n,1).^2;

nu = 0.1:0.1:50;

y = zeros(size(nu));
for i=1:length(nu)
  y(i) = sum(t2_lpdf(t2,nu(i)));
end

nu_ml = t_ml(t2,0.1)

figure
plot(nu,y)
  

function x = l1_rand(mu,sigma,M)
error('this function was just for testing')
nu = 1;
u = gamrnd(nu*1,2/nu,M);
mean_invu = mean(1./u)
x = normrnd(mu,sigma.*u,M);

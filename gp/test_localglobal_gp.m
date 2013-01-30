

N = 100;
Np = 10;

x = 1:N;               % data inputs
xp = linspace(1,N,Np); % pseudo inputs
xh = 1:(N+10);         % predictive inputs

logtheta_se = [log(1), log(N/10)];
logtheta_cs = log(3);

Kxpxp_se = gp_cov_se(xp,xp,logtheta_se
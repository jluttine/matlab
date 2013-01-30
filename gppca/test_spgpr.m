function test_spgpr

randn('state', 1)
rand('state', 2)

n = 1000;
m = 50;
d = 1;

s = 3;
inputs = 0.1*n*rand(d,n);
logtheta = -3;
logtheta = log(1.8411);
s2 = 0.03;
V = diag(sparse(ones(n,1))/s2);
y = colsum(sin(inputs)') + sqrt(s2)*randn(n,1);
c = V*y;

ind = randperm(n);
%pseudoinputs = 3*randn(d,m)+0;
pseudoinputs = sort(inputs(:,ind(1:m))); % random subset of inputs

LV = chol(V, 'lower');
[f, Covf, pseudoinputs_final, logtheta_final, Kpp] = gplearn([], V, inputs, ...
                                                  pseudoinputs, @gpcov, ...
                                                  logtheta, 'maxsearch', ...
                                                  10, 'cholv', LV, 'vy', c);

% $$$ figure
% $$$ plot3(inputs(1,:)', inputs(2,:)', V\c, 'ro')
% $$$ hold on
% $$$ plot3(pseudoinputs(1,:)', pseudoinputs(2,:)', ones(m,1), 'b+')
% $$$ plot3(pseudoinputs_final(1,:)', pseudoinputs(2,:)', -ones(m,1), 'k+')

%return

% $$$ figure
% $$$ pcolor(Covf);
% $$$ 
% $$$ figure
% $$$ pcolor(Kpp);

figure
inputsh = min(inputs-8):0.1:max(inputs+8);
[fh, varfh] = gppred(pseudoinputs_final, f, Covf, inputsh, @gpcov, ...
                     logtheta_final);
gpplot(inputsh, fh, 1*sqrt(varfh));
hold on

plot(inputs(1,:), V\c, 'ro')
plot(pseudoinputs(1,:), 1, 'b+')
plot(pseudoinputs_final(1,:), f, 'k+')

% $$$ theta = exp(logtheta)
% $$$ theta_final = exp(logtheta_final)

return;

return
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ N = 200;
% $$$ x = 1:N;
% $$$ 
% $$$ 
% $$$ % Covariance functions
% $$$ gpcov = @(in1,in2,theta) gpcov_ratquad(gpdist(in1,in2),theta(1), ...
% $$$                                        theta(2),theta(3));
% $$$ 
% $$$ % Covariance matrices
% $$$ theta = [1,20,1];
% $$$ Kxx = gpcov(x,x,theta);
% $$$ Kyy = Kxx + 0.1^2 * eye(N);
% $$$ 
% $$$ % Generate latent variables
% $$$ y = mvnrnd(zeros(1,N),Kyy)';
% $$$ %covfun = gpcov;
% $$$ %inittheta = theta;
% $$$ 
% $$$ % Noise level
% $$$ %s2 = 1 ^ 2;
% $$$ 
% $$$ %y = y + sqrt(s2)*randn(N,1);
% $$$ 
% $$$ xh = (min(x)-10):0.5:(max(x)+10);
% $$$ 
% $$$ %Kxx = gpcov(x,x,theta);
% $$$ %Kyy = Kxx + sq
% $$$ Kxhx = gpcov(xh,x,theta);
% $$$ Kxhxh = gpcov(xh,xh,theta);
% $$$ 
% $$$ % Full GP
% $$$ tic
% $$$ Kyy(Kyy < 0.1) = 0;
% $$$ R = chol(Kyy);
% $$$ yh = Kxhx * cholsolve(R, y);
% $$$ toc
% $$$ 
% $$$ figure
% $$$ plot(x, y, 'r+')
% $$$ hold on
% $$$ plot(xh, yh)
% $$$ 
% $$$ % Sparse GP
% $$$ tic
% $$$ Qyy = Kyy;
% $$$ %Coryy = Qyy;
% $$$ % Smoothly to zero?
% $$$ %ma = 0.12;
% $$$ %mi = 0.10;
% $$$ %Qyy(Qyy < ma & Qyy > mi) = (Qyy(Qyy < ma & Qyy > mi) - mi) * ;
% $$$ % $$$ Ryy = (1-Coryy) / (1-thres); % Scale thresold to 1
% $$$ % $$$ Ryy = Ryy ^ 5;
% $$$ % $$$ Ryy = Ryy / (1-thres);
% $$$ % $$$ Co = Coryy * 
% $$$ Qyy(Qyy < 0.1) = 0;
% $$$ Qyy = sparse(Qyy);
% $$$ [Rsp,p,S] = chol(Qyy);
% $$$ p
% $$$ yh = Kxhx * S * (Rsp \ (Rsp' \ (S'*y)));
% $$$ toc
% $$$ 
% $$$ figure
% $$$ plot(x, y, 'r+')
% $$$ hold on
% $$$ plot(xh, yh)

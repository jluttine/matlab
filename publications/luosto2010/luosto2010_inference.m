function luosto2010_inference

seed = 5
rand('state', seed);
randn('state', seed);

% $$$ N = 50;
% $$$ x = linspace(0.05, 0.95, N) + 0.05 * (rand(1,N)-0.5);
% $$$ %x = rand(1,N);
% $$$ 
% $$$ covfunc = @gpcov;
% $$$ real_logtheta = log(0.1);
% $$$ real_noise = 0.3;
% $$$ 
% $$$ % Generate data
% $$$ y = gprnd(x, real_logtheta, covfunc);
% $$$ yn = y + real_noise*randn(size(y));
% $$$ 
% $$$ 
% $$$ % Use Rasmussen's package for inference
% $$$ 
% $$$ % Normal
% $$$ ls = 0.1;
% $$$ mag = 1;
% $$$ noise = real_noise;
% $$$ logtheta = log([ls; mag; noise]);
% $$$ do_inference(x,yn,logtheta);
% $$$ 
% $$$ % Short
% $$$ ls = 0.03;
% $$$ mag = 1;
% $$$ noise = 0.0001;
% $$$ logtheta = log([ls; mag; noise]);
% $$$ do_inference(x,yn,logtheta);
% $$$ 
% $$$ % Long
% $$$ ls = 0.3;
% $$$ mag = 1;
% $$$ noise = 1;
% $$$ logtheta = log([ls; mag; noise]);
% $$$ do_inference(x,yn,logtheta);


N = 10;
x = linspace(0.05, 0.95, N) + 0.05 * (rand(1,N)-0.5);
%x = rand(1,N);

covfunc = @gpcov;
real_logtheta = log(0.1);
real_noise = 0.1;

% Generate data
y = gprnd(x, real_logtheta, covfunc);
% $$$ scale = sqrt(var(y,1));
% $$$ y = y / scale;
yn = y + real_noise*randn(size(y));


% Use Rasmussen's package for inference

% Normal
ls = 0.1;
mag = 1;
noise = real_noise;
logtheta = log([ls; mag; noise]);
do_inference(x,yn,logtheta, 'medium');

% Short
ls = 0.03;
mag = 1;
noise = 0.01;
logtheta = log([ls; mag; noise]);
do_inference(x,yn,logtheta, 'short');

% Long
ls = 0.3;
mag = 1;
noise = 0.5;
logtheta = log([ls; mag; noise]);
do_inference(x,yn,logtheta, 'long');


% $$$ % Normal
% $$$ ls = 0.1;
% $$$ mag = 1;
% $$$ noise = real_noise;
% $$$ logtheta = log([ls; mag; noise]);
% $$$ do_inference(x,yn,logtheta);
% $$$ 
% $$$ % Short
% $$$ ls = 0.03;
% $$$ mag = 1;
% $$$ noise = 0.01;
% $$$ logtheta = log([ls; mag; noise]);
% $$$ do_inference(x,yn,logtheta);
% $$$ 
% $$$ % Long
% $$$ ls = 0.3;
% $$$ mag = 1;
% $$$ noise = 0.5;
% $$$ logtheta = log([ls; mag; noise]);
% $$$ do_inference(x,yn,logtheta);
% $$$ 

function do_inference(x,yn, logtheta, name)
covfunc = {'covSum', {'covSEiso','covNoise'}};

% Optimize the parameters
xstar = 0:0.01:1;
%logtheta = minimize(logtheta, 'gpr', -100, covfunc, x(:), yn(:));
[mu S2] = gpr(logtheta, covfunc, x(:), yn(:), xstar(:));
S2 = S2 - exp(2*logtheta(3)); % minus the noise

figure
gpplot(xstar, mu, 2*sqrt(S2))
hold on
plot(x, yn, 'b+')

set_plot_style;
set(gca, 'xtick', [], 'ytick', []);
set_figure_size(gcf, 6, 5);

filename = sprintf('/home/jluttine/thesis/figures/fig_gpinference_%s', name);
print('-depsc2', filename);

% x = laplace_rand(mu, lambda, M, N)
%
% RANDL  random numbers distributed according to the Laplace distribution
%   RANDL(N) will return an NxN matrix containing pseudo-random values
%   drawn from a Laplace distribution with zero mean and standard deviation
%   one. RAND([M,N,...,P]) will return an MxNx...xP matrix.
%   RANDL([M,N,...,P],lambda) will return an MxNx...xP matrix of
%   pseudo-random numbers with parameter lambda. CAUTION: the pdf is
%   assumed as 
%           pdf = lambda/2 * exp(-lambda*abs(x-mu))
%
% The Laplace random numbers are generated using the the RAND function to
% generate uniformly distributed numbers in (0,1) and then the probability
% integral transformation, i.e. using the fact that 
%   if Fx(x) is the cdf of a random variable x, then the RV z=Fx(x) is
%   uniformly distributed in (0,1).
%
% In order to generate random numbers with mean mu and variance v, then
% generate the random numbers X with mean 0 and variance 1 and then
%       X = X*sqrt(v)+mu;
%
% C. Saragiotis, Oct 2008

function x = laplace_rand(mu, lambda, M, N)

if nargin < 1
  mu  = 0;
end
if nargin < 2 
  lambda = 1;
  % lambda = sqrt(2); % this gives a std=var=1.
end
if nargin < 3
  M = 1;
end
if nargin < 4
  N = 1;
end

z = rand(M,N);
x = zeros(M,N);

in = z<=.5;
ip = z> .5;
x(in) =  1/lambda *log(2*z(in));
x(ip) = -1/lambda *log(2*(1-z(ip)));


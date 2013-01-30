
% Test whether it is necessary to rotate VB PCA
function nc2010_rppca_toyexperiment(n,m,d,ncomps,rotate,maxiters)

if nargin < 3
  n = 200;
  m = 10;
  d = 4;
end

if nargin < 6
  maxiters = 10000;
end
if nargin < 5
  rotate = true;
end
if nargin < 4
  ncomps = d;
end

randn('state', 1);
rand('state', 1);

%error('Here is something wrong when generating the data!!')

% Generate multivariate normal data
mu = randn(m,1);
R = orth(randn(m,m)); % arbitrary rotation
eigs = (1 + [d:-1:1 zeros(1,m-d)]) .^ 2; % eigenvalues of covariance matrix
Cov = R*diag(eigs)*R';
Yn = mvnrnd(mu', Cov, n)';

% Divide into training and test sets
testset = (rand(m,n) < 0.2);
Ytrain = Yn;
Ytest = Yn;
Ytrain(testset) = NaN;
Ytest(~testset) = NaN;

% Generate some outliers to TRAINING set
p = (rand(1,n) < 0.02);
Ytrain(:,p) = 10 * randn(m,sum(p(:)));
Ytest(:,p) = nan; % remove outliers from test set
Ytrain(testset) = nan;


filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                    'nc2010_rppca_toyexperiment_n=%d_m=%d_d=%d_ncomps=' ...
                    '%d_rotate=%d'], n,m,d,ncomps, rotate);


results = rppcamv(Ytrain, ncomps, 'testset', Ytest, 'maxiters', maxiters, ...
                  'rotate', rotate, 'autosavetime', 1, 'autosavefile', ...
                  filename);
  

function nc2010_vbpca_mnistexperiment(flatW, rotate, ncomps)

M = load('mnist_all');

n = 100;

Y = double(M.train5');
inds = randperm(cols(Y));
Y = Y(:,inds(1:n));
size(Y)

m = rows(Y);

startupdatehyper = 1
maxiters = 1e4;

randn('state', 1);
rand('state', 1);

% Divide into training (90%) and test sets (10%)
testset = (rand(size(Y)) < 0.1);
Ytrain = Y;
Ytest = Y;
Ytrain(testset) = nan;
Ytest(~testset) = nan;

%size(Ytrain)

if flatW
  error('Using flat W might be problematic because m>n ???');
  disp('Using non-hierarchical flat prior for W');
  stringW = 'flatW';
else
  disp('Using hierarchical prior for W');
  stringW = 'hierW';
end

filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                    'nc2010_vbpca_mnistexperiment_%s_ncomps=%d_rotate=%d'], ...
                   stringW, ncomps, rotate);


results = vbpcamv(Ytrain, ncomps, 'testset', Ytest, 'maxiters', 1e7, ...
                  'rotate', rotate, 'startupdatehyper', startupdatehyper, ...
                  'startrotate', 1, 'autosavetime', 1, ...
                  'autosavefile', filename, 'fixw', flatW);
  

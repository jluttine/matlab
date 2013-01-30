
% Test whether it is necessary to rotate VB PCA
function nc2010_vbpca_mlexperiment(flatW, rotate, ncomps)

if nargin < 3
  ncomps = 100;
end

randn('state', 1);
rand('state', 1);

% Movie lens data with 100.000 ratings
U = load('/share/bayes/data/jaakko/movielens/u.data');
%Y = sparse(U(:,1), U(:,2), U(:,3));

% Divide into training (90%) and test sets (10%)
m = max(U(:,1));
n = max(U(:,2));
testset = (rand(rows(U),1) < 0.1);
Ytrain = full(sparse(U(~testset,1), U(~testset,2), U(~testset,3),m,n));
Ytest = full(sparse(U(testset,1), U(testset,2), U(testset,3),m,n));
Ytrain(Ytrain==0) = nan;
Ytest(Ytest==0) = nan;

%size(Ytrain)

if flatW
  disp('Using non-hierarchical flat prior for W');
  stringW = 'flatW';
else
  disp('Using hierarchical prior for W');
  stringW = 'hierW';
end

filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                    'nc2010_vbpca_mlexperiment_%s_ncomps=%d_rotate=%d'], ...
                   stringW, ncomps, rotate);


results = vbpcamv(Ytrain, ncomps, 'testset', Ytest, 'maxiters', 1e7, ...
                  'rotate', rotate, 'startupdatehyper', 1, ...
                  'startrotate', 1, 'autosavetime', 1, ...
                  'autosavefile', filename, 'fixw', flatW);
  

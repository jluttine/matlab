function nc2010_vbpca_compexperiment(flatW, datatype)

n = 200;
m = 50;
d = 10;

startupdatehyper = 2
maxiters = 1e5;

% Initialise random number generator
seed = mod(ceil(now),9999)
randn('state', seed);
rand('state', seed);

datasets = 10;

for ind=1:datasets

  fprintf('Dataset index = %d/%d\n', ind, datasets);
  
  % Generate multivariate normal data
  mu = randn(m,1);
  R = orth(randn(m,m)); % arbitrary rotation

  % eigenvalues of covariance matrix
  switch datatype
   case 1
    disp('Using weak subspace')
    eigs = (1 + [d:-1:1 zeros(1,m-d)]) .^ 2
    datastring = 'weak';
   case 2
    disp('Using strong subspace')
    eigs = ([5*ones(1,d) 1*ones(1,m-d)]) .^ 2
    datastring = 'strong';
   case 3
    disp('Using no subspace')
    eigs = (m:-1:1) .^ 2
    datastring = 'no';
  end  
  Cov = R*diag(eigs)*R';
  Yn = mvnrnd(mu', Cov, n)';

  % Generate some missing values
  Ytrain = Yn;
  Ytest = Yn;
  pmv = (rand(m,n) < 0.2);
  Ytrain(pmv) = NaN;
  Ytest(~pmv) = NaN;

  if flatW
    disp('Using non-hierarchical flat prior for W');
    stringW = 'flatW';
  else
    disp('Using hierarchical prior for W');
    stringW = 'hierW';
  end

  for ncomps = 1:m
    for rotate = [false true]
      filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                          'nc2010_vbpca_compexperiment_%s_subspace=%s_n=' ...
                          '%d_m=%d_d=%d_ncomps=%d_rotate=%d_seed=' ...
                          '%d_datasetindex=%d'], stringW, datastring, n,m, ...
                         d,ncomps, rotate, seed, ind);


      vbpcamv(Ytrain, ncomps, 'testset', Ytest, 'maxiters', maxiters, ...
              'rotate', rotate, 'startupdatehyper', startupdatehyper, ...
              'startrotate', 1, 'autosavetime', 1, 'autosavefile', filename, ...
              'fixw', flatW, 'convergence', 1e-6);
    end
  end
end
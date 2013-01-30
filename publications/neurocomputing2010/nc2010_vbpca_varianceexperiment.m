
% Test whether it is necessary to rotate VB PCA
function nc2010_vbpca_varianceexperiment(rotate, flatW, datatype)

debug = true;

n = 200;
m = 50;
d = 10;
ncomps = 30;

startupdatehyper = 1;
maxiters = 1000

randn('state', 1);
rand('state', 1);

% Generate multivariate normal data
if debug
  mu = zeros(m,1);
else
  mu = randn(m,1);
end
R = orth(randn(m,m)); % arbitrary rotation
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

% Generate some outliers (NO!)
Yno = Yn;

% Generate some missing values
missing_values = 0.2
Ynom = Yno;
pmv = (rand(m,n) < missing_values);
Ytest = Ynom;
Ynom(pmv) = NaN;
Ytest(~pmv) = NaN;

if flatW
  disp('Using non-hierarchical flat prior for W');
  stringW = 'flatW';
else
  disp('Using hierarchical prior for W');
  stringW = 'hierW';
end

% INITIALIZE LARGE/CORRECT NOISE VARIANCE
init.tau = 1e-2;

if false
  filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                      'nc2010_vbpca_varianceexperiment_%s_subspace=%s_n=' ...
                      '%d_m=%d_d=%d_ncomps=%d_rotate=%d'], stringW, ...
                     datastring,n,m,d,ncomps, rotate);
  autosavetime = 100;
else
  disp(['Debugging (using zero mean), saving the results to debug file - ' ...
        'press any key to continue!']);
  pause
  filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                      'nc2010_vbpca_variancedebugexperiment_%s_subspace=%s_n=' ...
                      '%d_m=%d_d=%d_ncomps=%d_rotate=%d'], stringW, ...
                     datastring,n,m,d,ncomps, rotate);
  mu = zeros(m,1);
  maxiters = 100;
  autosavetime = 20;
end


results = vbpcamv(Ynom, ncomps, 'testset', Ytest, 'maxiters', maxiters, ...
                  'rotate', rotate, 'startupdatehyper', startupdatehyper, ...
                  'startrotate', 1, 'autosavetime', autosavetime, 'autosavefile', ...
                  filename, 'fixw', flatW, 'init', init);
  

function [rmse_train, rmse_test] = aistats2012_mohsst5_rmse(F, testset, seed)

if nargin < 3
  [Ytrain, Ytest, Itrain, Itest, data, sea] = aistats2012_mohsst5_sets(testset);
else
  [Ytrain, Ytest, Itrain, Itest, data, sea] = aistats2012_mohsst5_sets(testset, ...
                                                    seed);
end

[M,N] = size(F);

% Error measures
[LON,LAT] = meshgrid(data.longitude, data.latitude);
weights = cosd(LAT(:));
weights = weights(sea);
weights = repmat(weights, [1, N]);
rmse_train = rmsew(Ytrain(Itrain)-F(Itrain),weights(Itrain));
rmse_test = rmsew(Ytest(Itest)-F(Itest),weights(Itest));
rmse_zero = rmsew(Ytest(Itest),weights(Itest));


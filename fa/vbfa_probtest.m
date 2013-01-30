% Compute the probability of test data for a trained VBFA model
% Usage:
%  [ prob, rms, Yrec, Ystd ] = vbfa_testprob( Q, Ytest )
%  [ prob, rms, Yrec, Ystd ] = vbfa_testprob( 'filename', Ytest )
function [ prob, rms, Yrec, Ystd ] = vbfa_probtest( Q, Ytest )

% Check that the test data coverage is different from the training
% data coverage
ObsTest = ~isnan(Ytest);

% Compute the reconstructions for the test data
Yrec = W'*X;

% Compute RMS

% Compute the probability of test data given the training data

% Compute the error bars

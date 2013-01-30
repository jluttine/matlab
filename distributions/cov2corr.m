% Transforms a covariance matrix into a correlation matrix

function C = cov2corr(C)

S = diag(diag(C).^(-0.5));
C = S*C*S;

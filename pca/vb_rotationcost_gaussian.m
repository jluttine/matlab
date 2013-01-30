
% P(X) = GAUSSIAN(MU0,COV0)
%
% Q(X) = GAUSSIAN(MU,COV)
%
% where the parameters are actually functions
%   MU = MU(THETA),
%   COV = COV(THETA).
%
% More generally, the whole joint distribution Q(X,MU0,COV0) can be a
% function of THETA.
%
% Evaluates KL(Q||P) w.r.t. THETA up to a constant, and the gradient.
%
% vb_rotationcost_gaussian(logdet_Cov, x_invCov0_x, x_invCov0_mu, ...
%                          mu_invCov0_mu, logdet_Cov0, grad_logdet_Cov, ...
%                          grad_x_invCov0_x, grad_x_invCov0_mu, ...
%                          grad_mu_invCov0_mu, grad_logdet_Cov0)
%
% LOGDET_COV         : <LOG(DET(COV))>
% X_INVCOV0_X        : <X'*INV(COV0)*X>
% X_INVCOV0_MU       : <X'*INV(COV0)*MU>
% MU_INVCOV0_MU      : <MU'*INV(COV0)*MU>
% LOGDET_COV0        : <LOG(DET(COV0))>
% GRAD_LOGDET_COV    : D(<LOG(DET(COV))>) / D(THETA)
% GRAD_X_INVCOV0_X   : D(<X'*INV(COV0)*X>) / D(THETA)
% GRAD_X_INVCOV0_MU  : D(<X'*INV(COV0)*MU>) / D(THETA)
% GRAD_MU_INVCOV0_MU : D(<MU'*INV(COV0)*MU>) / D(THETA)
% GRAD_LOGDET_COV0   : D(<LOG(DET(COV0))>) / D(THETA)
%
% where the expectations <...> are taken over the distribution
% Q(X,MU0,COV0). Note that the terms can be evaluated up to a constant
% w.r.t. THETA, which often gives huge savings.

% Last modified 2010-06-24
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function [KL,dKL] = vb_rotationcost_gaussian(logdet_Cov,         ...
                                             x_invCov0_x,        ...
                                             x_invCov0_mu,       ...
                                             mu_invCov0_mu,      ...
                                             logdet_Cov0,        ...
                                             grad_logdet_Cov,    ...
                                             grad_x_invCov0_x,   ...
                                             grad_x_invCov0_mu,  ...
                                             grad_mu_invCov0_mu, ...
                                             grad_logdet_Cov0)

% Cost from <log q(X)>
log_qX = -0.5 * logdet_Cov;
dlog_qX = -0.5 * grad_logdet_Cov;

% Cost from <log p(x)>
log_pX = -0.5 * (logdet_Cov0 + x_invCov0_x - 2*x_invCov0_mu + mu_invCov0_mu);
dlog_pX = -0.5 * (grad_logdet_Cov0 + grad_x_invCov0_x - 2*grad_x_invCov0_mu ...
                  + grad_mu_invCov0_mu);

KL = log_qX - log_pX;
dKL = dlog_qX - dlog_pX;




% $$$ function [f,df_dR] = vb_rotationcost_gaussian(R,logdet_R,invT_R, ...
% $$$                                               XX_R_invCov,X_Mu_invCov)
% $$$ 
% $$$ [D,N] = size(X);
% $$$ 
% $$$ % TODO: Check that the sizes match
% $$$ 
% $$$ % Cost from <log q(X)>
% $$$ log_qX = -N*logdet_R;
% $$$ dlog_qX = -N*invT_R;
% $$$ 
% $$$ % Cost from <log p(x)>
% $$$ log_pX = -0.5 * (traceprod(XX_R_invCov,R) - 2*traceprod(X_Mu_invCov,R));
% $$$ dlog_pX = X_Mu_invCov - XX_R_invCov;
% $$$ 
% $$$ f = log_pX -

%
% P(ALPHA) = GAMMA(A0,B0)
%
% Q(ALPHA) = GAMMA(A,B),
%
% where B = B(THETA).
%
% KL(Q||P) w.r.t. THETA up to a constant and the gradient.

function [KL,dKL] = vb_rotationcost_gamma(A0_logB,      ...
                                          A_B0_invB,    ...
                                          grad_A0_logB, ...
                                          grad_A_B0_invB)

KL = A0_logB + A_B0_invB;
dKL = grad_A0_logB + grad_A_B0_invB;
                                           
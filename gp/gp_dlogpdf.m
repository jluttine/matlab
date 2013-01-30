
% dl = gp_dlogpdf(y_invK_dK_invK_y, trace_invK_dK)

% dl = gp_dlogpdf(dK, invK, invK_y)

function dl = gp_dlogpdf(y_invK_dK_invK_y, trace_invK_dK)
dl = 0.5*y_invK_dK_invK_y - 0.5*trace_invK_dK;
end
% $$$ function dl = gp_dlogpdf(dK, invK, invK_y)
% $$$ dl = 0.5*(invK_y'*dK*invK_y) - 0.5*traceprod(invK,dK);
% $$$ end



% l = gp_logpdf(y_invK_y, logdetK, N)

function l = gp_logpdf(y_invK_y, logdetK, N)
l = gaussian_logpdf(y_invK_y, 0, 0, logdetK, N);
end

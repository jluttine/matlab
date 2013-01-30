%
% [K, dK] = gpK_decper(D, p1, p2, p3, p4)
%
% Decaying periodic covariance function for GP.
%
% p1 is the scale
% p2 is the decaying length
% p3 is the frequency of the periodicity
% p4 is the smoothness of the periodic signal
%
function [K, dK] = gpK_decper(D, p1, p2, p3, p4)
K = p1^2*exp(-0.5*(D.^2)/(p2^2)-2*sin(pi*D*p3).^2/(p4^2));

if nargout >= 2
  dK = zeros([size(D),4]);
  dK(:,:,1) = K .* 2 / p1;
  dK(:,:,2) = K .* (-0.5*D.^2) .* (-2*p2^(-3));
  dK(:,:,3) = K .* (-2*sin(pi*D*p3)*2) .* cos(pi*D*p3) .* (pi*D);
  dK(:,:,4) = K .* (-2*sin(pi*D*p3).^2) * (-2)*p4^(-3);
end

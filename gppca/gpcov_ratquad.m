
function [K,dK] = gpcov_ratquad(D, p1, p2, p3)
% f = 1 + D.^2/(2*p3^2*p2^2);
% K = p1^2*f.^(-p3^2);

f = 1 + D.^2/(2*p3^2*p2^2);
K = p1^2*f.^(-p3^2);

if nargout >= 2
  dK = zeros([size(D),3]);
  dK(:,:,1) = K .* 2 / p1;
  dK(:,:,2) = K .* (-p3^2).*f.^(-1) .* (-2).*D.^2/(2*p3^2)*p2^(-3);
  dK(:,:,3) = K .* (-2*p3*log(f) + (-p3^2)./f.*D.^2/(2*p2^2) * (-2) * p3^(-3));
end

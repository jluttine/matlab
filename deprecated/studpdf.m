function p = studpdf(x,mu,s2,nu)
warning('This function is deprecated')

% function p = studpdf(x,mu,s2,nu)

p = 1;
p = p .* gamma((nu+1)/2);
p = p ./ gamma(nu/2);
p = p ./ sqrt(nu.*pi.*s2);
p = p .* (1 + 1./nu .* ((x-mu).^2) ./ s2) .^ (-(nu+1)/2);

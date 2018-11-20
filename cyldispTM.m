function [detm, detd] = cyldispTM(cyl, epi)

% set epi if not provided
if nargin == 1
  epi = cyl.epi;
end

% interior exterior arguments
kperp = sqrt(cyl.k.^2.*cyl.ep.*cyl.mu - cyl.beta.^2);
kperpi = sqrt(cyl.k.^2.*epi.*cyl.mu - cyl.beta.^2);
ka = kperp.*cyl.a;
kia = kperpi.*cyl.a;

bessj = besselj(cyl.orders, kia);
bessjm = besselj(cyl.orders-1, kia);
bessjp = besselj(cyl.orders+1, kia);

% equation ratios
jrat =  0.5.*(bessjm-bessjp)./bessj.*kia;
hrat =  besselh1d(cyl.orders, ka)./besselh(cyl.orders, ka).*ka;

% derivatives
djratdk = (bessjp.*bessjm - bessj.*(bessjm-bessjp)./kia - bessj.^2)./kperpi./bessj.^2;
dkdepi = cyl.k.^2./2./kperpi;

% determinant derivative TODO rewrite
detd = jrat + epi.*djratdk.*dkdepi;

% determinant
detm = jrat-hrat;

newt = detm./detd;

function [detm,detd] = cyldispTE(cyl, epi)

% set epi if not provided
if nargin == 1
  epi = cyl.epi;
end

% interior exterior arguments
kperp = sqrt(cyl.k.^2.*cyl.ep.*cyl.mu ).*cyl.a;
kperpi = sqrt(cyl.k.^2.*epi.*cyl.mu ).*cyl.a;
ka = kperp;
kia = kperpi;

bessj = besselj(cyl.orders, kia);
bessjm = besselj(cyl.orders-1, kia);
bessjp = besselj(cyl.orders+1, kia);

% equation ratios
jrat = 0.5./kia.^2 + 0.5.*(bessjm-bessjp)./bessj./kia;
hrat = 0.5./ka.^2 + besselh1d(cyl.orders, ka)./besselh(cyl.orders, ka)./ka;

% derivatives %TODO rewrite
djratdk = (bessjp.*bessjm - bessj.*(bessjm-bessjp)./kia - bessj.^2)./kperpi./bessj.^2;
dkdepi = cyl.k.^2./2./kperpi;

% determinant derivative
detd = djratdk.*dkdepi;

 % determinant
detm = jrat-hrat;

%newt = detm./detd;

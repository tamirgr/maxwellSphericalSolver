function newt = sphrdispepinewt16(epi, sphr)

% evaluates ratio of dispersion relation and its derivative 
% for cylinder epsilon as variable

% set epi if not provided
if nargin == 1
  epi = sphr.epi;
end

% interior exterior arguments
% branch cut irrelevant
kperp = sqrt(sphr.a.^2.*sphr.k.^2.*sphr.ep);
kperpi = sqrt(sphr.a.^2.*sphr.k.^2.*epi);

ka = kperp;%maybe add .*sphr.a
kia = kperpi;%maybe add .*sphr.a

% scaled bessel functions
bessj = besselj(sphr.orders, kia, 1);
bessjm = besselj(sphr.orders-1, kia, 1);
bessjp = besselj(sphr.orders+1, kia, 1);

bessh = besselh(sphr.orders, 1, ka, 1);
besshm = besselh(sphr.orders-1, 1, ka, 1);
besshp = besselh(sphr.orders+1, 1, ka, 1);
bessj2m = besselj(sphr.orders-2, kia, 1);
bessj2p = besselj(sphr.orders+2, kia, 1);

% equation ratios
jrat = 1./kia.^2 + 0.5.*(bessjm-bessjp)./bessj./kia;
hrat = 1./ka.^2 +0.5.*(besshm - besshp)./bessh./ka;

% derivatives
djratdk = -1./kia/2 - 1./kia.^3 + ((bessjm - bessjp).*(-2*bessj + bessjp.*kia - bessjm.*kia) + kia.*(bessj2m.*bessj + bessj.*bessj2p))./4./(kia.*bessj).^2;
dkdepi = sphr.k.^2.*sphr.k.^2./kperpi/2;

detd = djratdk.*dkdepi;

% determinant
detm = (jrat-hrat);

newt = detm./detd;

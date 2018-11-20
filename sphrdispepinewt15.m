function newt = sphrdispepinewt15(epi, sphr)

% evaluates ratio of dispersion relation and its derivative 
% for cylinder epsilon as variable

% set epi if not provided
if nargin == 1
  epi = sphr.epi;
end

% interior exterior arguments
% branch cut irrelevant
kperp = sqrt(sphr.k.^2.*sphr.ep);
kperpi = sqrt(sphr.k.^2.*epi);

ka = kperp; %maybe add .*sphr.a
kia = kperpi; %maybe add .*sphr.a
orders = sphr.orders+0.5;
% scaled bessel functions
bessj = besselj(orders, kia, 1);
bessjm = besselj(orders-1, kia, 1);
bessjp = besselj(orders+1, kia, 1);
bessj2m = besselj(orders-2, kia, 1);
bessj2p = besselj(orders+2, kia, 1);

bessh = besselh(orders, 1, ka, 1);
besshm = besselh(orders-1, 1, ka, 1);
besshp = besselh(orders+1, 1, ka, 1);

% equation ratios
jrat =  0.5.*(bessjm-bessjp)./bessj.*kia;
hrat = 0.5.*(besshm-besshp)./bessh.*ka;

% derivatives
djratdk = -kia./2 + ((bessjm - bessjp).*(2*bessj + bessjp.*kia - bessjm.*kia) + kia.*(bessj2m.*bessj + bessj.*bessj2p))./4./(bessj.^2);
dkdepi = sphr.k.^2./kperpi/2;

% determinant derivative
detd = djratdk.*dkdepi;

% determinant
detm = (jrat-hrat);

newt = detm./detd;

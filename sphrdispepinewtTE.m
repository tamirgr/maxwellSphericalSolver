function newt = sphrdispepinewtTE(xnl, sphr)

% evaluates ratio of dispersion relation and its derivative 
% for cylinder epsilon as variable

% set epi if not provided
if nargin == 1
  xnl = sphr.epi;
end
l = sphr.orders;
% interior exterior arguments
% branch cut irrelevant
kperp = sphr.a.*sphr.k.*sqrt(sphr.ep);


% scaled bessel functions
bessj = SphericalBesselJ(l, xnl);
%bessjm = SphericalBesselJ(l-1, xnl);
bessjp = SphericalBesselJ(l+1, xnl);

% bessh = SphericalHenkelH(l, kperp);
% besshm = SphericalHenkelH(l-1, kperp);
% besshp = SphericalHenkelH(l+1, 1, kperp, 1);

constpart = DSphericalHankelH(l,kperp)/SphericalHankelH1(l, kperp)*kperp;


newt = (bessj.*((constpart-l).*bessj+xnl.*bessjp))./(xnl.*bessj.^2 - (2*l+1).*bessj.*bessjp+xnl.*bessjp.^2);
end
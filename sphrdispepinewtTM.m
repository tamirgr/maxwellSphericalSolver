function newt = sphrdispepinewtTM(xnl, sphr)

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

constpart = (kperp*DSphericalHankelH(l,kperp)/SphericalHankelH1(l, kperp)+1)/kperp^2;


newt = (xnl.*bessj.*(-(1+l - constpart*xnl.^2).*bessj + xnl.*bessjp))./ ...
    ((2+2*l+xnl.^2).*bessj.^2 - (2*l+3).*xnl.*bessj.*bessjp + xnl.^2.*bessjp.^2);


end

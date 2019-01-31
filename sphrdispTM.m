function [detm] = sphrdispTM(sphr, xnl)


if nargin == 1
  xnl = sphr.epi;
end

l = sphr.orders;
% interior exterior arguments
% branch cut irrelevant
kperp = sphr.a.*sphr.k.*sqrt(sphr.ep);


% equation ratios
jrat = (1 + DSphericalBesselJ(l,xnl).*xnl./SphericalBesselJ(l,xnl))./xnl.^2;
hrat = (kperp*DSphericalHankelH(l,kperp)/SphericalHankelH1(l, kperp)+1)/kperp^2;


 % determinant
detm = jrat-hrat;

%newt = detm./detd;
end
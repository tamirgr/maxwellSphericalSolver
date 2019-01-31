function [detm] = sphrdispTE(sphr, xnl)


if nargin == 1
  xnl = sphr.epi;
end

l = sphr.orders;
% interior exterior arguments
% branch cut irrelevant
kperp = sphr.a.*sphr.k.*sqrt(sphr.ep);


% equation ratios
jrat = DSphericalBesselJ(l,xnl).*xnl./SphericalBesselJ(l,xnl);
hrat = DSphericalHankelH(l,kperp)/SphericalHankelH1(l, kperp)*kperp;



 % determinant
detm = jrat-hrat;

%newt = detm./detd;

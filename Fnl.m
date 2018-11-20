%r is a vector, returns a vector
function [value] = Fnl(r,sphr,epiNL,n,l)
%'a' is the radius of the sphere
ext = r < sphr.a;
tmp = ext.*r + sphr.a.*(~ext);

xnl = (sphr.k.*sqrt(epiNL(l+1,n)./sphr.ep)).*tmp;
% A = normalizationCoeffCalc(sphr,xnl,mode);

% A = sqrt(sphr.a./r).*(~ext)./besselh(l+0.5,sphr.k*sphr.a).*besselh(l+0.5,sphr.k.*r) + ext;
A = 1; % ONLY FOR DEBUG!!!
value = SphericalBesselJ(l,xnl).*A;

end
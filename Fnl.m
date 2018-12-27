function [value] = Fnl(r,sphr,epiNL,n,l)
% TO DO - function should receive the mode (TE or TM)
ext = r < sphr.a;
tmp = r.*ext + sphr.a.*(~ext);

xnl = (sphr.k.*sqrt(epiNL(l+1,n)./sphr.ep)).*r;
% A = normalizationCoeffCalc(sphr,xnl,'M');
A = normalizationCoeffCalc(sphr,xnl,'E'); %This line configures configured for TE field. 

% A = sqrt(sphr.a./r).*(~ext)./besselh(l+0.5,sphr.k*sphr.a).*besselh(l+0.5,sphr.k.*r) + ext;
% A = 1; % ONLY FOR DEBUG!!!
value = SphericalBesselJ(l,xnl).*A;

end
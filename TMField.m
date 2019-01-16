function [R,Th,Phi] = TMField(r,th,phi,sphr,epiNL,n,l,m)
%% Input Check
if l-1<m
    R   = 0;
    Th  = 0;
    Phi = 0;
	return
end
%% Init
ext = (r<=sphr.a);
rho = sphr.k * sqrt(epiNL(l,n)/sphr.ep) .* r;
xnl = sphr.k * sphr.a * sqrt(epiNL(l,n)/sphr.ep); % 
coeffJH = SphericalBesselJ(l-1,xnl)/SphericalHankelH1(l-1,sphr.k*sphr.a);
normce = normalizationCoeffCalc(l-1,sphr,rho,'M');
%% Calculation
[resRIn,resThIn,resPhiIn] = curlFXlmPoint2(r,th,phi,epiNL,sphr,n,l-1,m);
[resROut,resThOut,resPhiOut] = curlFXlmPointHenkel2(r,th,phi,epiNL,sphr,n,l-1,m);    
%% Result
R   = 1i/sphr.k*(ext.*resRIn/epiNL(l,n)*sphr.ep   + coeffJH*(~ext).*resROut).*normce;
Th  = 1i/sphr.k*(ext.*resThIn/epiNL(l,n)*sphr.ep  + coeffJH*(~ext).*resThOut).*normce;
Phi = 1i/sphr.k*(ext.*resPhiIn/epiNL(l,n)*sphr.ep + coeffJH*(~ext).*resPhiOut).*normce;
end
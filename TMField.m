function [R,Th,Phi] = TMField(r,th,phi,sphr,epiNL,n,l,m)
%% Input Check
if l<m
    R   = 0;
    Th  = 0;
    Phi = 0;
	return
end
%% Init
int = (r<=sphr.a);
if size(epiNL)>1
    epi = epiNL(l+1,n);
else
    epi = epiNL;
end
rho = sphr.k * sqrt(epi/sphr.ep) .* r;
xnl = sphr.k * sphr.a * sqrt(epi/sphr.ep); % 
coeffJH = SphericalBesselJ(l,xnl)/SphericalHankelH1(l,sphr.k*sphr.a);
normce = normalizationCoeffCalc(l,sphr,rho,'M');
%% Calculation
[resRIn,resThIn,resPhiIn] = curlFXlmPoint2(r,th,phi,epiNL,sphr,n,l,m);
[resROut,resThOut,resPhiOut] = curlFXlmPointHenkel2(r,th,phi,epiNL,sphr,n,l,m);    
%% Result
R   = 1i/sphr.k*(int.*resRIn/epi*sphr.ep   + coeffJH*(~int).*resROut).*normce;
Th  = 1i/sphr.k*(int.*resThIn/epi*sphr.ep  + coeffJH*(~int).*resThOut).*normce;
Phi = 1i/sphr.k*(int.*resPhiIn/epi*sphr.ep + coeffJH*(~int).*resPhiOut).*normce;
end
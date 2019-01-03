function [R,Th,Phi] = TMField(r,th,phi,sphr,epiNL,n,l,m)
%% Input Check
if l<m
    R   = 0;
    Th  = 0;
    Phi = 0;
	return
end
%% Init
ext = (r<=sphr.a);

xnl = sphr.k * sphr.a * sqrt(epiNL(l,n)/sphr.ep); % 
coeffJH = SphericalBesselJ(l,xnl)/SphericalHankelH1(l,sphr.k*sphr.a);
%% Calculation
[resRIn,resThIn,resPhiIn] = curlFXlmPoint2(r,th,phi,epiNL,sphr,n,l,m);
[resROut,resThOut,resPhiOut] = curlFXlmPointHenkel2(r,th,phi,epiNL,sphr,n,l,m);    
%% Result
R   = 1i/sphr.k*(ext.*resRIn/epiNL(l,n)*sphr.ep   + coeffJH*(~ext).*resROut);
Th  = 1i/sphr.k*(ext.*resThIn/epiNL(l,n)*sphr.ep  + coeffJH*(~ext).*resThOut);
Phi = 1i/sphr.k*(ext.*resPhiIn/epiNL(l,n)*sphr.ep + coeffJH*(~ext).*resPhiOut);
end
function [norm] = normalizationCoeffCalc(l,sphr,xnl,mode)
% normalizationCoeffCalc calculates the value of field coeffiecients A or
% B.

if mode == 'M'
    normA = 1;
    elA = sphr.a.^3/(2.*xnl.^2).*(...
        psitag(l,xnl).^2 + psi(l,xnl).^2 - l.*(l+1).*SphericalBesselJ(l,xnl).^2 + psitag(l,xnl).*SphericalBesselJ(l,xnl));
elseif mode == 'E'
%     normA = xnl.^3./sphr.a.^3;
    normA = 1;
%     elA = 1/2.*xnl.^3.*(SphericalBesselJ(l,xnl).^2 - SphericalBesselJ(l-1,xnl).*SphericalBesselJ(l+1,xnl));
    elA = sphr.a.^3/(2.*xnl.^2).*(...
        psitag(l,xnl).^2 + psi(l,xnl).^2 - l.*(l+1).*SphericalBesselJ(l,xnl).^2 - psitag(l,xnl).*SphericalBesselJ(l,xnl));
else
    error('element mode in SphereGeometry is used incorrectly. check normalizationCoeffCalc.m');
end

% norm = 1;
norm = sqrt(normA./elA);
% norm = (elA+elB+elC); %this line is just to debug the value of the integrals.
end

function [res] =  psi(n,p)
    res = p.*SphericalBesselJ(n, p);
end

function [res] =  psitag(n,p)
    res = SphericalBesselJ(n, p) + ...
    p.*(-(SphericalBesselJ(n, p)./(2*p)) + 0.5*(SphericalBesselJ(-1 + n, p) - SphericalBesselJ(1 + n, p)));
end
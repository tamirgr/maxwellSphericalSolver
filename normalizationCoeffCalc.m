function [norm] = normalizationCoeffCalc(sphr,xnl,mode)
% normalizationCoeffCalc calculates the value of field coeffiecients A or
% B.
l = sphr.orders;

if mode == 'E'
    normA = 2*xnl./sphr.a.^3/pi;
    elA = SphericalBesselJ(l,xnl).*DSphericalBesselJ(l,xnl).*xnl.^2;
    elB = 1/4*pi.*xnl.^2.*(besselj(l+0.5,xnl).^2 - besselj(l-0.5,xnl).*besselj(l+1.5,xnl));
    elC = (-l*(l+1) + ((l+1)*l+1)/((4*l+2)*xnl.^2) ) .* calcSquaredBesselIntegral(l,xnl); 
%    elC = (-l*(l+1) ) .* calcSquaredBesselIntegral(l,xnl); 
elseif mode == 'M'
    normA = xnl.^3./sphr.a.^3;
    elA = 0;
%     elB = 1/4*pi.*xnl.^2.*(besselj(l+0.5,xnl).^2 - besselj(l-0.5,xnl).*besselj(l+1.5,xnl));
    elB = 1/2.*xnl.^3.*(SphericalBesselJ(l,xnl).^2 - SphericalBesselJ(l-1,xnl).*SphericalBesselJ(l+1,xnl));
    elC = 0;
else
    error('element mode in SphereGeometry is used incorrectly. check normalizationCoeffCalc.m');
end

norm = sqrt(normA./(elA+elB+elC));
% norm = (elA+elB+elC); %this line is just to debug the value of the integrals.
end

function [res] = calcSquaredBesselIntegral(l, a)
coeff = pi*a.^(2*l+1)/(2*l+1)^3/(4^l*gamma(l+0.5)^2);

res = coeff * hypergeom([l+0.5,l+1],[l+1.5,l+1.5,2*l+2],-a.^2);
end
function [norm] = normalizationCoeffCalc(l,sphr,xnl,mode)
% normalizationCoeffCalc calculates the value of field coeffiecients A or
% B.

if mode == 'E'
    normA = 2*xnl./sphr.a.^3/pi;
    elA = SphericalBesselJ(l,xnl).*DSphericalBesselJ(l,xnl).*xnl.^2;
    elB = pi/4.*xnl.^2.*(SphericalBesselJ(l,xnl).^2 - SphericalBesselJ(l-1,xnl).*SphericalBesselJ(l+1,xnl));
    elC = (-l.*(l+1) + ((l+1).*l+1)./((4*l+2).*xnl.^2) ) .* calcSquaredBesselIntegral(l,xnl); 
%    elC = (-l*(l+1) ) .* calcSquaredBesselIntegral(l,xnl); 
elseif mode == 'M'
    normA = xnl.^3./sphr.a.^3;
    elA = 0;
%     elB = 1/4*pi.*xnl.^2.*(besselj(l+0.5,xnl).^2 - besselj(l-0.5,xnl).*besselj(l+1.5,xnl));
    elB = 1/2.*xnl.^3.*(SphericalBesselJ(l,xnl).^2 - SphericalBesselJ(l-1,xnl).*SphericalBesselJ(l+1,xnl));
    % perhaps this should be divided by 4 instead of 2
    elC = 0;
else
    error('element mode in SphereGeometry is used incorrectly. check normalizationCoeffCalc.m');
end

norm = sqrt(normA./(elA+elB+elC));
% norm = (elA+elB+elC); %this line is just to debug the value of the integrals.
end

function [res] = calcSquaredBesselIntegral(l, a)

for i = 1:length(a)
    for j = 1:length(a)
        for k = 1:length(a)
            v= linspace(0,a(i,j,k),100);
            res = trapz(v,forIntegral(v,l));
        end
    end
end
%coeff = pi*a.^(2*l+1)/4^(l+1)*factorial(2*l);
%res = coeff .* 1;%hypergeom([l+0.5,l+1],[l+1.5,l+1.5,2*l+2],-a.^2);
end
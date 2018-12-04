function [ res ] = DSphericalBesselJ( l , r)
%DSphericalBesselJ calculate the first derivative of a spherical bessel
%function of the first kind at point r of order l.
res = -1./(2.*r).*SphericalBesselJ(l,r) + 0.5*(SphericalBesselJ(l-1,r) - SphericalBesselJ(l+1,r) );

end


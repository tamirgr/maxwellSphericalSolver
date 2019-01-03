function [ res ] = DSphericalHankelH( l , r)
%DSphericalBesselJ calculate the first derivative of a spherical bessel
%function of the third kind at point r of order l.
res = -1./(2.*r).*SphericalHankelH1(l,r) + 0.5*(SphericalHankelH1(l-1,r) - SphericalHankelH1(l+1,r) );

end


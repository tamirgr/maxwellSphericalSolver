function [ res ] = SphericalBesselJ( nu, z)
%SphericalBesselJ is the spherical bessl function of the first kind
% nu - order of bessel function. nu must be a positive integer.
% z - value at which to evaluate.

res = sqrt(pi/2./z).*besselj(nu+0.5, z);
end


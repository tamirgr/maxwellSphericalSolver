function [ res ] = SphericalHankelH1( nu, z)
%SphericalHankelH1 is the spherical bessl function of the third kind
% nu - order of bessel function. nu must be a positive integer.
% z - value at which to evaluate.

res = sqrt(pi/2./z).*besselh(nu+0.5, z);
end
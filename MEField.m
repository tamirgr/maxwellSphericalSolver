function [resR, resTh,resPhi] = MEField(r,th,phi,sphr,epiNL,n,l,m)
rho = sphr.k * sqrt(epiNL(l,n)/sphr.ep) .* r;

coeff = 1/sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial

resR = zeros(size(r));

resTh = -m*sin(m.*phi).*coeff./plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*SphericalBesselJ(l,rho).*normalizationCoeffCalc(sphr,rho,'E');

resPhi = -cos(m.*phi).*coeff./plmcoeff.*(...
    m.*cot(th).*legendrePlm(l,m,cos(th))...
    +legendrePlm(l,m+1,cos(th))...
    ).* SphericalBesselJ(l,rho).*normalizationCoeffCalc(sphr,rho,'E');
end
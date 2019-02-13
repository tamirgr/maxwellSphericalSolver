function [resR,resTh,resPhi] = NOField(r,th,phi,sphr,l,m)
%rho = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;
rho = sphr.k * r;
normalizationcoeff = normalizationCoeffCalc(l,sphr,rho,'E');
coeff = 1/sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = 1;%sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial

resR = l*(l+1)*sin(m*phi).*plmcoeff.*legendrePlm(l,m,cos(th)).*SphericalBesselJ(l,rho)./rho.*normalizationcoeff;

resPhi = cos(m*phi)*m.*coeff.*plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*(SphericalBesselJ(l,rho)./rho + DSphericalBesselJ(l,rho)).*normalizationcoeff;

resTh = sin(m*phi)*coeff.*plmcoeff.*(...
    m.*cot(th).*legendrePlm(l,m,cos(th))...
    +legendrePlm(l,m+1,cos(th))...
    ).* (SphericalBesselJ(l,rho)./rho + DSphericalBesselJ(l,rho)).*normalizationcoeff;
end
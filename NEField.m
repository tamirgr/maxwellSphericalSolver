function [resR,resTh,resPhi] = NEField(r,th,phi,sphr,epiNL,n,l,m)
rho = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;
normalizationcoeff = normalizationCoeffCalc(l,sphr,rho,'M');
coeff = 1/sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = 1;%sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial

resR = l*(l+1)*cos(m*phi).*plmcoeff.*legendrePlm(l,m,cos(th)).*SphericalBesselJ(l,rho)./rho.*normalizationcoeff;
 
resPhi = -m*sin(m*phi).*coeff.*plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*(SphericalBesselJ(l,rho)./rho + DSphericalBesselJ(l,rho)).*normalizationcoeff;

resTh = cos(m*phi).*coeff.*plmcoeff.*(...
    m.*cot(th).*legendrePlm(l,m,cos(th))...
    +legendrePlm(l,m+1,cos(th))...
    ).* (SphericalBesselJ(l,rho)./rho + DSphericalBesselJ(l,rho)).*normalizationcoeff;
end
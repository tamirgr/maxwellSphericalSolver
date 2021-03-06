function [resR, resTh,resPhi] = MOField(r,th,phi,sphr,mu,l,m)
%rho = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;

rho = sphr.k *mu *r;
normalizationcoeff = normalizationCoeffCalc(l,sphr,rho,'E');
coeff = 1/sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = 1;%sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial

resR = zeros(size(r));

resTh = cos(m*phi).* m.*coeff.*plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*SphericalBesselJ(l,rho).*normalizationcoeff;

resPhi = -sin(m*phi).*coeff.*plmcoeff.*(...
    m.*cot(th).*legendrePlm(l,m,cos(th))...
    +legendrePlm(l,m+1,cos(th))...
    ).* SphericalBesselJ(l,rho).*normalizationcoeff;
end
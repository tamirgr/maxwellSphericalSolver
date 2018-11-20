function [xlmTh,xlmPhi] = TEField(r,th,phi,sphr,epiNL,n,l,m)
rho = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;

coeff = exp(1i.*m.*phi)./sin(th).*sqrt(pi/(4*l+2)); % result of calculation with normalization coefficient for spherical bessel function
ylmcoeff = sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
ylmcoeff1 = sqrt(4*pi/(2*l+3)*factorial(l+m+1)/factorial(l-m+1)); % normalization coefficient for associated legendre ploynomial

xlmTh = 1i.*coeff.*ylmcoeff.*legendrePlm(l,m,cos(th)).*SphericalBesselJ(l,rho);

xlmPhi = coeff.*(...
    (l+1).*cos(th).*ylmcoeff.*legendrePlm(l,m,cos(th)) +...
    (l-m+1).*ylmcoeff1.*legendrePlm(l+1,m,cos(th))...
    ).* SphericalBesselJ(l,rho);
end
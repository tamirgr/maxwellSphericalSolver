function [xlmTh,xlmPhi] = TEField(r,th,phi,sphr,epiNL,n,l,m)
rho = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;

coeff = 1i*exp(1i.*m.*phi)./sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial

xlmTh = 1i.*m.*coeff./plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*SphericalBesselJ(l,rho);

xlmPhi = coeff./plmcoeff.*(...
    legendrePlm(l,m+1,cos(th))...
    +m.*cot(th).*legendrePlm(l,m,cos(th))...
    ).* SphericalBesselJ(l,rho);
end
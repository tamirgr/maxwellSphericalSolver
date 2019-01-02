function [xlmR,xlmTh,xlmPhi] = TMField(r,th,phi,sphr,epiNL,n,l,m)
rho = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;

coeff = 1i*exp(1i.*m.*phi)./sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial

xlmTh = 1i.*m.*coeff./plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*normalizationCoeffCalc(sphr,rho,'E');

xlmPhi = -1*coeff./plmcoeff.*(...
    m.*cot(th).*legendrePlm(l,m,cos(th))...
    +legendrePlm(l,m+1,cos(th))...
    ).*normalizationCoeffCalc(sphr,rho,'E');

xnl = sphr.k * sphr.a * sqrt(epiNL(l+1,n)/sphr.ep); % 
coeffJH = SphericalBesselJ(l,xnl)/SphericalHankelH1(l,sphr.k*sphr.a);

xlmTh  = xlmTh.*SphericalBesselJ(l,rho).*(r<=sphr.a)  + coeffJH*xlmTh.*SphericalHankelH1(l,rho).*(r>sphr.a);
xlmPhi = xlmPhi.*SphericalBesselJ(l,rho).*(r<=sphr.a) + coeffJH*xlmPhi.*SphericalHankelH1(l,rho).*(r>sphr.a);
end
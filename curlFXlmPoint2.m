function [rv, thetav, phiv] = curlFXlmPoint2(r, th, phi, epiNL, sphr,n,l,m)
% curlFXlmPoint2 calculates the rotor of the multiplication of f_nl times
% X_lm in a secific coorediante given by r, th and phi.
% alpha is the permittivity coefficient of the length r.

%% Input Check
if l<m
    rv = 0;
	thetav = 0;
    phiv = 0;
    return
end

if size(epiNL)>1
    epi = epiNL(l+1,n);
else
    epi = epiNL;
end

%% Init
if nargin < 5
    %sphr = SphereGeometry;
    r = 1.5;
    th = pi/6;
    phi = pi/4;
    l = 3;
    m=2;
    n=4;
    %epiNL = 0.5;
    %sphr.k = 3.0;
    %sphr.ep = 4.0;
    %alphar = sphr.k * (sphr.ep - epiNL) / sphr.epi * r;
    rho = 3.0 * (4.0 - 0.5) / 4.0 * 1.5;
else
%    alphar = sphr.k * (sphr.ep - epiNL) / sphr.ep * r;

%     A = sqrt(sphr.a./r).*(~ext)./besselh(l+0.5,sphr.k*sphr.a).*besselh(l+0.5,sphr.k.*r) + ext;
    rho = sphr.k * sqrt(epi/sphr.ep) .* r;
end

%% Coefficients

ylmcoeff = sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m));

coeff = exp(1i.*m.*phi)./rho;
coeffang =  coeff.*csc(th);

coefftheta =  -1*coeffang;
coeffphi =  1i.*m.*coeffang;

coeffsb = (l+1) .* SphericalBesselJ(l,rho) - rho .* SphericalBesselJ(l+1,rho);
%% calculation
 
rv = coeff .* SphericalBesselJ(l,rho) * l * (l+1) ./ ylmcoeff .* legendrePlm(l,m,cos(th));

thetav  = coefftheta .* (...
    m .* cos(th) ./ ylmcoeff .* legendrePlm(l,m,cos(th)) +  sin(th)./ ylmcoeff .* legendrePlm(l,m+1,cos(th))...
    ) .* coeffsb;

phiv= coeffphi ./ ylmcoeff .* legendrePlm(l,m,cos(th)) .* coeffsb;
    
end
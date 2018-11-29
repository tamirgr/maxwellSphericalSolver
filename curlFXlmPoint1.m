function [rv, thetav, phiv] = curlFXlmPoint1(r, th, phi, epiNL, sphr,n,l,m)
% curlFXlmPoint calculates the rotor of the multiplication of f_nl times
% X_lm in a secific coorediante given by r, th and phi.
% alpha is the permittivity coefficient of the length r.

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
    kunlR = 3.0 * (4.0 - 0.5) / 4.0 * 1.5;
else
%    alphar = sphr.k * (sphr.ep - epiNL) / sphr.ep * r;

%     A = sqrt(sphr.a./r).*(~ext)./besselh(l+0.5,sphr.k*sphr.a).*besselh(l+0.5,sphr.k.*r) + ext;
    kunlR = sphr.k * sqrt(epiNL(l+1,n)/sphr.ep) .* r;
end

%% Coefficients

dcoeff =  sqrt(2) /( sqrt(l*(l+1)))./kunlR;
ylm = Ylm(l,m,th,phi);

ylmP1 = Ylm(l,m+1,th,phi);
ylmP2 = Ylm(l,m+2,th,phi);

coeffsb =  (- kunlR .* SphericalBesselJ(l+1,kunlR) + (l+1).*SphericalBesselJ(l,kunlR));

%these 3 coefficients are calculated here to make the formula more readable.
coeffr1 = m*(m+1);
coeffr2 = 2 * (m+1) * sqrt((l-m)*(l+m+1)) .* cot(th) .* exp(-1i * phi);
coeffr3 = exp(-2*1i .* phi) .* sqrt((l-m-1)*(l-m)*(l+m+1)*(l+m+2));
 
%% calculation
 
rv = dcoeff .* SphericalBesselJ(l,kunlR) .* ...
    ( coeffr1 .* ylm - coeffr2 .* ylmP1 - coeffr3 .* ylmP2 );

thetav = -1.* dcoeff .* coeffsb .* ...
    (  m .* cot(th) .* ylm + sqrt((l-m)*(l+m+1)) .*  exp(-1i * phi)  .* ylmP1 );

phiv = dcoeff .* (1i) .* m .* csc(th) .* ylm .* coeffsb;
end
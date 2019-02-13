function [R, Th,Phi] = TEField(r,th,phi,sphr,epiNL,n,l,m)
R=0;
%% Input Check
if l<m
    Th = 0;
	Phi = 0;
    return
end
%% Init
if size(epiNL)>1
    epi = epiNL(l+1,n);
else
    epi = epiNL;
end
rho = sphr.k * sqrt(epi/sphr.ep) .* r;
ext = (r<=sphr.a);
coeff = 1i*exp(1i.*m.*phi)./sqrt(l*(l+1)); % result of calculation with normalization coefficient for spherical bessel function
plmcoeff = 1.0;%sqrt(4*pi/(2*l+1)*factorial(l+m)/factorial(l-m)); % normalization coefficient for associated legendre ploynomial
% plmcoeff1 = sqrt(4*pi/(2*l+1)*factorial(l+m+1)/factorial(l-m-1)); % normalization coefficient for associated legendre ploynomial
normce = normalizationCoeffCalc(l,sphr,rho,'E');
%% Calculation
Th = 1i.*m.*coeff./plmcoeff./sin(th).*legendrePlm(l,m,cos(th)).*normce;

Phi = -1*coeff./plmcoeff.*(...
    m.*cot(th).*legendrePlm(l,m,cos(th))...
    +legendrePlm(l,m+1,cos(th))...
    ).*normce;

xnl = sphr.k * sphr.a * sqrt(epi/sphr.ep); % 
coeffJH = SphericalBesselJ(l,xnl)/SphericalHankelH1(l,sphr.k*sphr.a);

spbj = SphericalBesselJ(l,rho);
sph1 = SphericalHankelH1(l,rho);

%% Result
Th  = Th.*spbj.*ext  + coeffJH*Th.*sph1.*(~ext);
Phi = Phi.*spbj.*ext + coeffJH*Phi.*sph1.*(~ext);
end
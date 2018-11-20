function [rv, thetav, phiv] = curlFXlmPointHenkel(r, th, phi, epiNL, sphr)
% curlFXlmPoint calculates the rotor of the multiplication of f_nl times
% X_lm in a secific coorediante given by r, th and phi.
% alpha is the permittivity coefficient of the length r.
if nargin < 5
    %sphr = SphereGeometry;
    r = 1.5;
    th = pi/6;
    phi = pi/4;
    l = 3;
    %epiNL = 0.5;
    %sphr.k = 3.0;
    %sphr.ep = 4.0;
    %alphar = sphr.k * (sphr.ep - epiNL) / sphr.epi * r;
    alphar = 3.0 * (4.0 - 0.5) / 4.0 * 1.5;
else
    alphar = sphr.k *   sqrt(epiNL/sphr.ep) * r;
end
alpha = sphr.k *   sqrt(epiNL/sphr.ep) * sphr.a;

rv = zeros([size(epiNL),2*size(epiNL,1)+1]);
phiv = zeros([size(epiNL),2*size(epiNL,1)+1]);
thetav = zeros([size(epiNL),2*size(epiNL,1)+1]);


for l = 0:(size(epiNL,1)-1) 
    m = -l:l;

    dcoeff = 1 /(2 * sqrt(l*(l+1)))/r;
    ylm = cat(3,Ylm(l,th,phi),[0],[0]);


    coeffang = besselj(l,alpha(l+1,:))./besselh(l,sphr.a*sphr.k).* sphr.k .* ( besselh(l-1,alphar(l+1,:)) - besselh(l+1,alphar(l+1,:))) ...
                    + 2 .* besselh(l,alphar(l+1,:));

    max(coeffang)
    %% calculation
    for mm = 1:(2*l+1)
        %these 3 coefficients are calculated here to make the formula more
        %readable.
        coeffr1 = m(mm)*(m(mm)+1);
        coeffr2 = 2 * cot(th) * exp(-1i * phi) * ((m(mm)+1) * sqrt((l-m(mm)) *(l+m(mm)+1)));
        coeffr3 = exp(-2i * phi) * sqrt((l-m(mm))*(l-m(mm)+1)*(l+m(mm)+1)*(l+m(mm)+2));

        rv(l+1,:,mm) = besselj(l,alpha(l+1,:))./besselh(l,sphr.a*sphr.k).* dcoeff .* besselh(m(mm),alphar(l+1,:)) .* ...
            ( coeffr1 .* ylm(1,1,mm) - ...
            coeffr2 .* ylm(1,1,mm+1) - ...
            coeffr3 .* ylm(1,1,mm+2) );
        phiv(l+1,:,mm) = 1i .*  exp(-1i * phi) .* dcoeff .* coeffang .* ...
            ( exp(1i .* phi) .* m(mm) .* cot(th) .* ylm(1,1,mm) + ...
            sqrt((l-m(mm))*(l+m(mm)+1)) .* ylm(1,1,mm+1) );
   %     phiv(1,:,1)
        thetav(l+1,:,mm) = dcoeff .* (-1) .* m(mm) .* csc(th) .* ylm(1,1,mm) .* coeffang;
    end
end




end
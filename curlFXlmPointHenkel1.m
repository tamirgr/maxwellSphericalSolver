function [rv, thetav, phiv] = curlFXlmPointHenkel1(r, th, phi, epiNL, sphr,n,l,m)
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
    alphar = sphr.k *   sqrt(epiNL(l+1,n)/sphr.ep) * r;
end
alpha = sphr.k *   sqrt(epiNL(l+1,n)/sphr.ep) * sphr.a;

%for l = 0:(size(epiNL,1)-1) 
   % m = -l:l;

    dcoeff = 1 /(2 * sqrt(l*(l+1)))./r;
    ylm = Ylm(l,m,th,phi);
    if(m+1<=l)
        ylmP = Ylm(l,m+1,th,phi);
        if(m+2<=l)
            ylmP2 = Ylm(l,m+2,th,phi);
        else
            ylmP2 = zeros(size(th));
        end
    else
        ylmP = zeros(size(th));
         ylmP2 = zeros(size(th));
    end
    coeffang = besselj(l,alpha)./besselh(l,sphr.a*sphr.k).* sphr.k .* ( besselh(l-1,alphar) - besselh(l+1,alphar)) ...
                    + 2 .* besselh(l,alphar);

   % max(coeffang)
    %% calculation
  %  for mm = 1:(2*l+1)
        %these 3 coefficients are calculated here to make the formula more
        %readable.
        coeffr1 = m*(m+1);
        coeffr2 = 2 * cot(th) .* exp(-1i * phi) * ((m+1) * sqrt((l-m) *(l+m+1)));
        coeffr3 = exp(-2i * phi) * sqrt((l-m)*(l-m+1)*(l+m+1)*(l+m+2));

        rv = besselj(l,alpha)./besselh(l,sphr.a*sphr.k).* dcoeff .* besselh(m,alphar) .* ...
            ( coeffr1 .* ylm - ...
            coeffr2 .* ylmP - ...
            coeffr3 .* ylmP2 );
        phiv = 1i .*  exp(-1i * phi) .* dcoeff .* coeffang .* ...
            ( exp(1i .* phi) .* m .* cot(th) .* ylm + ...
            sqrt((l-m)*(l+m+1)) .* ylmP );
   %     phiv(1,:,1)
        thetav = dcoeff .* (-1) .* m .* csc(th) .* ylm .* coeffang;
        
   % end
%end




end
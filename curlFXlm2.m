function [rv, thetav, phiv] = curlFXlm2(r, th, phi, l, alpha)
% curlFXlm2 calculates the rotor of the multiplication of f_nl times X_lm
% alpha is the permittivity coefficient of the length r.
% this is the testing function for curlFXlm2 which calculates each point
% seperately, rather than multiplying matrices.
tic
if nargin < 5
    lr = 50;
    r = linspace(0,10,lr);
    lth = 40;
    th = linspace(0,pi,lth);
    lphi = 30;
    phi = linspace(-pi,pi,lphi);
    l = 3;
    alpha = 0.5;
else
    lr = length(r);
    lth = length(th);
    lphi = length(phi);
end
%% init
m = -l:l;
alphar = alpha.*r;
tmpcoeff = 1 /(2 * sqrt(l*(l+1)));
dcoeff = ones(lr, lth, lphi, length(m)) * tmpcoeff;
ylm = Ylm(l,th,phi);

rv = zeros(lr, lth, lphi, length(m));
phiv = zeros(lr, lth, lphi, length(m));
thetav = zeros(lr, lth, lphi, length(m));

for mm = 1:length(m)
    coeffr1 =  m(mm)*(m(mm)+1);
    for rr = 1:lr
        dcoeff(rr,:,:,:) = bsxfun(@rdivide, dcoeff(rr,:,:,:), r(rr) ); % used to simplify result of curl since it's a joined component
        coeffang = alphar(rr) * (besselj(l-1,alphar(rr)) - besselj(l+1,alphar(rr))) ...
                    + 2 * besselj(l,alphar(rr));
        for tt = 1:lth
            for pp = 1:lphi
                
                tmprv = 0.0+1i*0;
                tmpphiv = 0.0+1i*0;
                if (m(mm)+1)<=l
                    coeffr2 = 2 * cot(th(tt)) * exp(-1i * phi(pp)) * ((m(mm)+1) * sqrt((l-m(mm)) *(l+m(mm)+1)));
                    tmprv = coeffr2 * ylm(tt,pp,mm+1);
                    tmpphiv = ylm(tt,pp,mm+1);
                    if (mm+2)<=l
                        coeffr3 = exp(-2i * phi(pp)) * sqrt((l-m(mm))*(l-m(mm)+1)*(l+m(mm)+1)*(l+m(mm)+2));
                        tmprv = tmprv + coeffr3 * ylm(tt,pp,mm+2);
                    end
                end
                rv(rr,tt,pp,mm) = dcoeff(rr,tt,pp,mm) * besselj(m(mm),alphar(rr)) * ...
                    ( coeffr1 * ylm(tt,pp,mm) - tmprv ); 
                phiv(rr,tt,pp,mm) = 1i *  exp(-1i * phi(pp)) * dcoeff(rr,tt,pp,mm) * coeffang * ...
                    ( exp(1i * phi(pp)) * m(mm) * cot(th(tt)) * ylm(tt,pp,mm) + ...
                    sqrt((l-m(mm))*(l+m(mm)+1)) *  tmpphiv) ;
                thetav(rr,tt,pp,mm) = dcoeff(rr,tt,pp,mm) * -m(mm) * csc(th(tt)) *...
                    ylm(tt,pp,mm) * coeffang;
            end
        end
    end
end

toc
end
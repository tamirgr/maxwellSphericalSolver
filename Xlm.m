function [xlmTh,xlmPhi]= Xlm(th,phi,l,m)

xlmPhi = -1i/sqrt(l*(l+1)).*dYlm(l,m,th,phi);

xlmTh = -1*m*Ylm(l,m,th,phi)./sin(th)/sqrt(l*(l+1));
end
function [xlmTh,xlmPhi]= Xlm(th,phi,l,m)

xlmTh = -1*m*Ylm(l,m,th,phi)./sin(th)/sqrt(l*(l+1));

xlmPhi = (-1i)/sqrt(l*(l+1)).*dYlm(l,m,th,phi);
end
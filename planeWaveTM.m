function [ res ] = planeWaveTM( H0, l, a ,coeff1, coeff2)
%planeWaveTM 
%
fac = 2*pi*1i^(l+1)*2*(l+1)*l; %is this the real coefficient in the 
res = fac*integrateJ2R2(l ,a ,coeff1, coeff2);

res = H0*res;

end


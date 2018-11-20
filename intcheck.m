function [delta] = intcheck(a,l)
    x=linspace(a/1000,a,1000);
    I1 = trapz(x,forIntegral(x,l))
    I2 = pi/2*a^2/2*(besselj(l,a,1)^2-besselj(l+1,a,1)*besselj(l-1,a,1));
    %I2 = pi/2*a^2/2*(besselj(l,a,1)^2-besselj(l+1,a,1)*besselj(l-1,a,1))
    
    delta = abs(I1-I2)/I2;
end
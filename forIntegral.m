function [v] = forIntegral(x,l)
 v = besselj(l,x).*besselj(l,x).*x*pi/2;
 
 %v =pi/2*(besseljd(l,x).^2.*x - besselj(l,x,1).*besseljd(l,x) + besselj(l,x,1).^2./x/4);  
end
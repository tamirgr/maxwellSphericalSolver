function [ res ] = integrateJ2R2( l ,a ,coeff1, coeff2)
%integrateJ2R2 
% l - bessel index number
% a - radius of sphere (for integration)
% coeff1 - length coefficient (k*sqrt(epsilon))
% coeff2 - second length coefficient for the second bessel function
fac = a.*pi./2./sqrt(coeff1.*coeff2)/(coeff1.^2-coeff2.^2);
res = fac.*(coeff2.*besselj(-(1/2) + l, a.*coeff2).* besselj(3/2 + l, a.*coeff1) - ...
            coeff1.*besselj(-(1/2) + l, a.*coeff1).* besselj(3/2 + l, a.*coeff2));
end
% if l1 == l2
% else
%     res = 2.^(-3 - l1 - l2).* a.^(3 + l1 + l2).*pi.* gamma(2 + l1 + l2).*...
%         gamma(1/2 *(3 + l1 + l2)) .* hypergeom([1/2 *(2 + l1 + l2),1/2 *(3 + l1 + l2), 1/2 *(3 + l1 + l2)],...
%             [3/2 + l1, 1/2 *(5 + l1 + l2), 3/2 + l2, 2 + l1 + l2],-a.^2);
% end


function [value] = dYlm(l,m, th, phi)

ylm = Ylm(l,m, th, phi);

x = cot(th);
%m = -l:l;
%value = zeros(size(th));

if(m~=l)
    ylmP = Ylm(l,m+1,th,phi);
    value = m*x.*ylm + ...
        sqrt((l-m+1)*(l+m+2))*exp(-1i*phi).*ylmP;
else
    value = l*x.*ylm; 
end
end
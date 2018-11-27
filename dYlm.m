function [value] = dYlm(l,m, th, phi)

ylm = Ylm(l,m, th, phi);
ylmP = Ylm(l,m+1,th,phi);

value = m.*cot(th).*ylm + sqrt((l-m).*(l+m+1)) .* exp(-1i*phi).*ylmP;
end
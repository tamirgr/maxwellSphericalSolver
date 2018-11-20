function [normA] = normalizationCoeffCalc(sphr,xnl,modes)
normA = sqrt(2/sphr.a^3)
l = sphr.orders+0.5;
if(modes == 'M')
    normA = normA./besselj(l+1,xml,1);
else
    normA = normA.*sqrt(xnl.*shpr.*sphr.epi)./sqrt(pi*(besselj(l,xml,1).^2 - besselj(l-1,xml,1).*besselj(l+1,xml,1))+numericalIntegration(sphr)); 
    %normA = normA.*sqrt(xnl.*shpr.*sphr.epi)./sqrt(pi*(besselj(l,xml,1).^2 - besselj(l-1,xml,1).*besselj(l+1,xml,1)+(l^2 + 0.75)/(2*l)*(besselj(l,xml,1).^2./xnl-(djl(l-1,xlm).*besselj(l,xnl,1) - djl(l,xlm).*besselj(l-1,xnl,1))./xnl.^2)));
end
end
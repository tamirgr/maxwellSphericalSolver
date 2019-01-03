function [a,b,c,d]=calcCoeff(sphr,epiNL,n,l)
    m = sqrt(epiNL(1:l,n));
    x = sphr.k.*sphr.a*ones(l,1);
    meu0 = sphr.ep;
    meu1 = sqrt(1.5);
    tmp1 = meu0*m.^2; 
    a = (tmp1.*pJJ(m.*x,x,transpose(1:l)) - meu1*pJJ(x,m.*x,transpose(1:l)))./(tmp1.*pJH(x,m.*x,transpose(1:l)) - meu1*pHJ(x,m.*x,transpose(1:l))); 
    b = (meu1.*pJJ(m.*x,x,transpose(1:l)) - meu0.*pJJ(x,m.*x,transpose(1:l)))./(meu1.*pJH(x,m.*x,transpose(1:l)) - meu0.*pHJ(x,m.*x,transpose(1:l))); 
    
    c = (meu1.*pJH(x,x,transpose(1:l)) - meu1.*pHJ(x,x,transpose(1:l)))./(meu1.*pJH(x,m.*x,transpose(1:l)) - meu0.*pHJ(x,m.*x,transpose(1:l))); 
    d = (meu1.*m.*pJH(x,x,transpose(1:l)) - meu1.*m.*pHJ(x,x,transpose(1:l)))./(meu0.*m.^2.*pJH(x,m.*x,transpose(1:l)) - meu1.*pHJ(x,m.*x,transpose(1:l))); 

end

function [res] = pJH(x,y,l)
     res = SphericalBesselJ(l,y).*(SphericalHankelH1(l,x) + x.*DSphericalHankelH(l,x));
end 

function [res] = pHJ(x,y,l)
     res = SphericalHankelH1(l,x).*(SphericalBesselJ(l,y) + y.*DSphericalBesselJ(l,y));
end

function [res] = pJJ(x,y,l)
     res = SphericalBesselJ(l,x).*(SphericalBesselJ(l,y) + y.*DSphericalBesselJ(l,y));
end
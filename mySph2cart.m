function [X,Y,Z] = mySph2cart(rField,thField,phiField,th,phi)
    s1= sin(th);
    s2 = sin(phi);
    c1 = cos(th);
    c2 = cos(phi);
    X = s1.*c2.*rField + c1.*c2.*thField - s2.*phiField;
    Y = s1.*s2.*rField + c1.*s2.*thField + c2.*phiField;
    Z = c1.*rField - s1.*thField ;
end

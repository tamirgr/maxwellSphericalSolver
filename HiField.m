function [resR,resTh,resPhi] = HiField(r,th,phi,sphr,epiNL,n,l,m,HiR,HiTh,HiPhi,H0r,H0th,H0phi,coeff1,coeff2)
    [Mr,Mth,Mphi] = MEField(r,th,phi,sphr,epiNL,n,l,m);
    [Nr,Nth,Nphi] = NOField(r,th,phi,sphr,epiNL,n,l,m);

    resR = HiR   - sqrt(epiNL(l,n)).*sphr.k/sphr.ep.*H0r.*1i^l.*sqrt((2.*l + 1)/l/(l + 1)).*(coeff2(l).*Mr + 1i.*coeff1(l).*Nr);
    resTh = HiTh  - sqrt(epiNL(l,n)).*sphr.k/sphr.ep.*H0th.*1i^l.*sqrt((2.*l + 1)/l/(l + 1)).*(coeff2(l).*Mth + 1i.*coeff1(l).*Nth);
    resPhi = HiPhi  - sqrt(epiNL(l,n)).*sphr.k/sphr.ep.*H0phi.*1i^l.*sqrt((2.*l + 1)/l/(l + 1)).*(coeff2(l).*Mphi + 1i.*coeff1(l).*Nphi);
end
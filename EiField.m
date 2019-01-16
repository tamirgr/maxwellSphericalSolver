function [resR,resTh,resPhi] = EiField(l,EiR,EiTh,EiPhi,E0r,E0th,E0phi,coeff1,coeff2,Mr,Mth,Mphi,Nr,Nth,Nphi)

    resR = EiR + E0r.*1i^l.*sqrt((2.*l + 1)/l/(l + 1)).*(coeff1(l).*Mr - 1i.*coeff2(l).*Nr);
    resTh = EiTh + E0th.*1i^l.*sqrt((2.*l + 1)/l/(l + 1)).*(Mth.*coeff1(l) - 1i.*coeff2(l).*Nth);
    resPhi = EiPhi + E0phi.*1i^l.*sqrt((2.*l + 1)/l/(l + 1)).*(Mphi.*coeff1(l) - 1i.*coeff2(l).*Nphi);
end
function value = Ylm(l, m, th, phi)
if l>= abs(m)
    norm = sqrt((2.*l+1)./4./pi./factorial(l+m).*factorial(l-m));
    elev = legendrePlm(l, m, cos(th));
    % azim = cos(m.*phi);
    azim = exp(1i.*m.*phi);

    value = norm .* elev .* azim;
else
    value = 0;
end
end

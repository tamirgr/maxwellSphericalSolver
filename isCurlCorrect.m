sphr =SphereGeometry;

% background properties
sphr.ep = 1; %sphr.mu = 1.0;

% field properties
sphr.k = 1.0; %sphr.beta = 0.0;

% sphrinder properties
sphr.a = 1.0;
%sphr.mui = 1;
%sphr.mu = 2;

% sphrinder coordinates
len = 60;
range = 2;
sphr.x = 0.0; sphr.y = 0.0;
sphr.z = 0.0;
x = linspace(-range,range,len);
y = linspace(-range,range,len);
z = linspace(-range,range,len);
l = 1; %l indx1
n = 1; %n indx2
m = 1; %m indx3
sphr.ordersN = 20;
%sphr.beta = 0.5;

% epiNL1 = zeros(l,nmodes);
epiNL2 = zeros(l,sphr.ordersN);
for k = 1:10
    sphr.orders = k;
%     epiNL1(k,:) = disprootsepi2(sphr, nmodes);
    epiNL2(k,:) = disprootsepi3(sphr, sphr.ordersN);
end
sphr.orders = l; 
[X,Y,Z] = meshgrid(x,y,z);
[phi,th,r] = cart2sph(X,Y,Z);
th = pi/2 - th;
U = X./X*0.1;
fnl = Fnl(r, sphr, epiNL2, n, l);
[xlmTh,xlmPhi] = Xlm(th,phi,l, m);
% coeffIn = 1i/sphr.k/(1-sphr.ep);
% coeffOut = 1i/sphr.k;
[cx, cy, cz] = mySph2cart(0,fnl.*xlmTh,fnl.*xlmPhi,th,phi);
% [cx, cy, cz] = sph2cart(0,fnl.*xlmTh,fnl.*xlmPhi);
% [eRCurl, eThCurl, ePhiCurl]= curl(fnl,xlmTh,xlmPhi);
% [ex, ey, ez, cav] = curl(X,Y,Z,cx, cy, cz);
[ex, ey, ez, cav] = curl(cx, cy, cz);
kz = 1;
ex = ex.*kz;
ey = ey.*kz;
ez = ez.*kz;
% [eRCurlreal, eThCurlreal, ePhiCurlreal] = cart2sph(real(ex), real(ey), real(ez));
% [eRCurlimg, eThCurlimg, ePhiCurlimg] = cart2sph(imag(ex), imag(ey), imag(ez));
% eRCurl = eRCurlreal + 1i*eRCurlimg;
% eThCurl = eThCurlreal + 1i*eThCurlimg;
% ePhiCurl = ePhiCurlreal + 1i*ePhiCurlimg;
[resRIn1,resThIn1,resPhiIn1] = curlFXlmPoint1(r,th,phi,epiNL2,sphr,l,n,m);
[resRIn2,resThIn2,resPhiIn2] = curlFXlmPoint2(r,th,phi,epiNL2,sphr,l,n,m);
[resX1, resY1, resZ1] = mySph2cart(resRIn1,resThIn1,resPhiIn1,th,phi);
[resX2, resY2, resZ2] = mySph2cart(resRIn2,resThIn2,resPhiIn2,th,phi);
% [resX, resY, resZ] = sph2cart(resRIn,resThIn,resPhiIn);
% [resROut,resThOut,resPhiOut] = curlFXlmPointHenkel1(r,th,phi,epiNL2,sphr,l,n,m);    
% rsub =   resRIn -  eRCurl;
% phisub =   resPhiIn -  ePhiCurl;
% thsub =   resThIn -  eThCurl;
% xsub = (ex - real(resX))./real(resX);
% ysub = (ey - real(resY))./real(resY);
% zsub = (ez - real(resZ))./real(resZ);
xsub1 = (ex - resX1)./resX1;
ysub1 = (ey - resY1)./resY1;
zsub1 = (ez - resZ1)./resZ1;
xsub2 = (ex - resX2)./resX2;
ysub2 = (ey - resY2)./resY2;
zsub2 = (ez - resZ2)./resZ2;
xsub3 = (resX1 - resX2)./resX2;
ysub3 = (resY1 - resY2)./resY2;
zsub3 = (resZ1 - resZ2)./resZ2;
% xsub4 = (ex - resX2)./max(max(max(resX2)));
% ysub4 = (ey - resY2)./max(max(max(resY2)));
% zsub4 = (ez - resZ2)./max(max(max(resZ2)));
% xsub = (ex - real(resX))./ex;
% ysub = (ey - real(resY))./ey;
% zsub = (ez - real(resZ))./ez;
% subres1 = [max(max(max(xsub1))), min(min(min(xsub1))) ;
% max(max(max(ysub1))), min(min(min(ysub1))) ;
% max(max(max(zsub1))), min(min(min(zsub1)))]
% subres2 = [max(max(max(xsub2))), min(min(min(xsub2))) ;
% max(max(max(ysub2))), min(min(min(ysub2))) ;
% max(max(max(zsub2))), min(min(min(zsub2)))]
subres3 = [max(max(max(xsub3))), min(min(min(xsub3))) ;
max(max(max(ysub3))), min(min(min(ysub3))) ;
max(max(max(zsub3))), min(min(min(zsub3)))]
% subres4 = [max(max(max(xsub4))), min(min(min(xsub4))) ;
% max(max(max(ysub4))), min(min(min(ysub4))) ;
% max(max(max(zsub4))), min(min(min(zsub4)))]
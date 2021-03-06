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
l = 2; %l indx1
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

fnl = Fnl(r, sphr, epiNL2, n, l,'E');
[xlmTh,xlmPhi] = Xlm(th,phi,l, m);
xlmTh1 = xlmTh.*fnl;
xlmPhi1 = xlmPhi.*fnl;

[xlmTh2,xlmPhi2] = TEField(r,th,phi,sphr,epiNL2,n,l,m);

subphi1 = (xlmPhi1 - xlmPhi2)./xlmPhi2;
subth1 = (xlmTh1 - xlmTh2)./xlmTh2;
subres1 = [max(max(max(subphi1))), min(min(min(subphi1))) ;
max(max(max(subth1))), min(min(min(subth1)))]

[ExRot,EyRot,EzRot] = mySph2cart(0,xlmTh2,xlmPhi2,th,phi);
    
ExR = real(ExRot);
EyR = real(EyRot);
EzR = real(EzRot);
     
dispx = [x(floor(len/3)),x(floor(len/3*2))];
    
figure;
slice(X,Y,Z,ExR,dispx,dispx,dispx);
colorbar();
shading interp
title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));


% figure;
% hist(real(xlmTh1(:)),800);
% title('xlmTh.*fnl');
% figure;
% hist(real(xlmTh2(:)),800);
% title('TEField theta');
% figure;
% hist(real(xlmPhi1(:)),800);
% title('xlmPhi.*fnl');
% figure;
% hist(real(xlmPhi2(:)),800);
% title('TEField phi');
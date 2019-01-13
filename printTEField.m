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
len = 40;
range = 2;
sphr.x = 0.0; sphr.y = 0.0;
sphr.z = 0.0;
x = linspace(-range,range,len);
y = linspace(-range,range,len);
z = linspace(-range,range,len);
l = 2; %l indx1
n = 10; %n indx2
m = 2; %m indx3
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

% [R,Th,Phi] = TEField(r,th,phi,sphr,epiNL2,n,l,m);
[Jx, Jy, Jz] = genWave(len, range, 'x plane polarized', 0.5, range, sphr.k);
epsinc = sqrt(1.5);
epsback = 1.0;
ExZero = zeros(len,len,len);
EyZero = zeros(len,len,len);
EzZero = zeros(len,len,len);
[R,Th,Phi] = genTEField( ExZero, EyZero, EzZero, Jx, Jy, Jz, epsinc, epsback, sphr, n, l, m);

[ExRot,EyRot,EzRot] = mySph2cart(R,Th,Phi,th,phi);
    
ExR = real(ExRot);
EyR = real(EyRot);
EzR = real(EzRot);
     
dispx = [x(floor(len/3)),x(floor(len/3*2))];
    
figure;
colormap('jet');
slice(X,Y,Z,ExR,dispx,dispx,dispx);
colorbar();
shading interp
% caxis([-0.004    0.004]);
title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));

figure;
colormap('jet');
slice(X,Y,Z,EyR,dispx,dispx,dispx);
colorbar();
shading interp
% caxis([-0.004    0.004]);
title(sprintf('EyReal n=%d, l=%d, m=%d',n,l,m));

figure;
colormap('jet');
slice(X,Y,Z,EzR,dispx,dispx,dispx);
colorbar();
shading interp
% caxis([-0.004    0.004]);
title(sprintf('EzReal n=%d, l=%d, m=%d',n,l,m));
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
l = 3; %l indx1 relative position in epiNL. subtract 1 for non-epiNL use.
n = 3; %n indx2
m = 0; %m indx3
sphr.ordersN = 20;
%sphr.beta = 0.5;

sphr.orders = l-1; 
[X,Y,Z] = meshgrid(x,y,z);
[phi,th,r] = cart2sph(X,Y,Z);
th = pi/2 - th;

% [R,Th,Phi] = TMField(r,th,phi,sphr,epiNL2,n,l,m);
[Jx, Jy, Jz] = genWave(len, range, 'x plane polarized', 1.0, 0, sphr.k);
epsinc = sqrt(1.5);
epsback = 1.0;
ExZero = 1i./sphr.k.*1.*Jx;
EyZero = 1i./sphr.k.*1.*Jy;
EzZero = 1i./sphr.k.*1.*Jz;
% ExZero = Jx;
% EyZero = Jy;
% EzZero = Jz;
[ExRot,EyRot,EzRot] = genTMField( ExZero, EyZero, EzZero, Jx, Jy, Jz, epsinc, epsback, sphr, sphr.ordersN, l, m);
    
ExR = real(ExRot);
EyR = real(EyRot);
EzR = real(EzRot);
     
% dispx = [x(floor(len/3)),x(floor(len/3*2))];
dispx = [x(floor(len/2))];
    
figure;
colormap('jet');
slice(X,Y,Z,ExR,dispx,dispx,dispx);
colorbar();
shading interp
% caxis([-abs(max(max(max(ExR))))    abs(max(max(max(ExR))))]);
% caxis([-1    1]*1*10^-6);
title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));
% alpha(0.5);

% figure;
% colormap('jet');
% slice(X,Y,Z,EyR,dispx,dispx,dispx);
% colorbar();
% shading interp
% % caxis([-abs(max(max(max(EyR))))    abs(max(max(max(EyR))))]);
% title(sprintf('EyReal n=%d, l=%d, m=%d',n,l,m));
% 
% figure;
% colormap('jet');
% slice(X,Y,Z,EzR,dispx,dispx,dispx);
% colorbar();
% shading interp
% % caxis([-abs(max(max(max(EzR))))    abs(max(max(max(EzR))))]);
% title(sprintf('EzReal n=%d, l=%d, m=%d',n,l,m));
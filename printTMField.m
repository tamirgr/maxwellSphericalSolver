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
% [Jx, Jy, Jz] = genWave(len, range, 'x plane polarized', 1.0, 0, sphr.k);
% epsinc = sqrt(1.5);
epsback = 1.0;

[R,Th,Phi] = genTMField2( epsback, sphr, sphr.ordersN, l, m, len, range);
[ExRot,EyRot,EzRot] = mySph2cart(R,Th,Phi,th,phi);

ExR = real(ExRot);
EyR = real(EyRot);
EzR = real(EzRot);
     
% dispx = [x(floor(len/3)),x(floor(len/3*2))];
dispx = [x(floor(len/2))];

displayFields( ExR , EyR , EzR ,X,Y,Z, n,l,m,dispx);    
% 
% figure;
% colormap('jet');
% slice(X,Y,Z,ExR,dispx,dispx,dispx);
% colorbar();
% shading interp
% % caxis([-abs(max(max(max(ExR))))    abs(max(max(max(ExR))))]);
% % caxis([-1    1]*1*10^-6);
% title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));
% % alpha(0.5);
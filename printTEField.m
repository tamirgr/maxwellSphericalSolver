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
n = 20; %n indx2
m = 1; %m indx3
sphr.ordersN = 20;
%sphr.beta = 0.5;
mu = ones(l+1,n)*(sqrt(1.5));

sphr.orders = l; 
[X,Y,Z] = meshgrid(x,y,z);
[phi,th,r] = cart2sph(X,Y,Z);
th = pi/2 - th;

% [R,Th,Phi] = TEField(r,th,phi,sphr,epiNL2,n,l,m);
% [Jx, Jy, Jz] = genWave(len, range, 'x plane polarized', 1.0, 0.0, sphr.k);
% epsinc = sqrt(1.5);
epsback = 1.0;
close all
[ExRot,EyRot,EzRot] = genTEField2(epsback, sphr, n, l, m, len, range);
   
ExR = real(ExRot);
EyR = real(EyRot);
EzR = real(EzRot);

dispx = [x(floor(len/2))];
    
displayFields( ExR , EyR , EzR ,X,Y,Z, n,l,m,dispx,1, 0);

[Ex,Ey,Ez,Hx,Hy,Hz] = MieElectricalFields(sphr, x, y, z,l,n,m,ones(size(r)),ones(size(r)),ones(size(r)),0,0,0,mu);
displayFields( real(Ex) , real(Ey) , real(Ez) ,X,Y,Z, n,l,m,dispx,3, 0);

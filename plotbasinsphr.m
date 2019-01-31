function []=plotbasinsphr(order,k,r)
% plots example attraction basins and contour selection
factor = k;
sphr = SphereGeometry;

% background properties
sphr.ep = 1; %sphr.mu = 1;

% field properties
sphr.k = k; sphr.a = r;

% cylinder properties
%sphr.epi = 1; sphr.mui = 1;
sphr.orders = order;
 
%sphr.beta = 0.0; sphr.phi = 0;

a = -20; b = 25;
reals = linspace(a, b, 1000);
c = -15; d = 15;
imags = linspace(c, d, 500);

[re, im] = meshgrid(reals, imags);

epis = re + 1i.*im;

detm = sphrdispTM(sphr, epis);

% plot attraction basins
figure(2)
[~, chan] = contourf(re, im, log10(abs(detm)));
caxis([-2,2]);
% adjust level curves and color map
set(chan, 'LevelList', -2.5:0.2:6.5);
end
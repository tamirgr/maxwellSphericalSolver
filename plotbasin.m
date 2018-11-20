function []=plotbasin(order,k)
% plots example attraction basins and contour selection
factor = k;
cyl = CylinderGeometry;

% background properties
cyl.ep = 1; cyl.mu = 1;

% field properties
cyl.k = factor; cyl.a = 1;

% cylinder properties
cyl.epi = 2; cyl.mui = 1;
cyl.orders = order;

cyl.beta = 0.0; cyl.phi = 0;

a = -10/factor; b = 250/factor;
reals = linspace(a, b, 1000);
c = -20/factor; d = 15/factor;
imags = linspace(c, d, 500);

[re, im] = meshgrid(reals, imags);

epis = re + i.*im;

detm = cyldispTM(cyl, epis);

% plot attraction basins
figure(2)
[~, chan] = contourf(re, im, log10(abs(detm)));

% adjust level curves and color map
set(chan, 'LevelList', [-2.5:0.2:6.5])
end
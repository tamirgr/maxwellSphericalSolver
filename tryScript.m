function []=tryScript(c)
% define a sphrinder
sphr =SphereGeometry;

% background properties
sphr.ep = 1; sphr.mu = 1.0;

% field properties
sphr.k = 1.0; sphr.beta = 0.0;

% sphrinder properties
sphr.a = 1.0;
sphr.mui = 1;
sphr.orders = 0.0;
fintol = 1e-5;
%roots = apm(@sphrdispepinewt15, 120-1.6*i, 5, 1, 0, [], [], fintol, sphr);
t = 0:0.01:2*pi;
r=1;
%c = 22-1.6*i;
expt = exp(1i*t);
  z = c + r*expt;
temp = 1./sphrdispepinewt16(z,sphr).*expt;
res = r/length(t)*sum(temp)
proot = newton(@sphrdispepinewt16, c, 5, sphr)

zero = abs(sphrdispepinewt16(proot,sphr))

end
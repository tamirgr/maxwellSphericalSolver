function [roots, epis, cen, rad] = disprootsepi3(sphr, nroots)

% Finds roots of the sphrinder dispersion relation using argument principle method 
%
% Usage
% roots = disprootsepi(sphr, nroots)
% [roots, sings, cen, rad] = disprootsepi(sphr, nroots)
%
% Inputs
% sphr - structure containing sphrinder parameters
% nroots - number of roots desired, beginning from fundamental
%
% Outputs
% roots - locations of roots
% sings - (optional) locations of singularities of dispersion relation
% cen - (optional) centers of contours
% rad - (optianal) radii of contours

% number of contours
% each singularity known to have 2 roots nearby
ncntrs = ceil(nroots);

% get zero of J, corresponds to kperpi*a
jzero = besselzero(sphr.orders+0.5, max(4, ncntrs+1));
% calculate corresponding epsilon
%epis = ((jzero./sphr.a).^2 + sphr.beta.^2)./sphr.k.^2./sphr.mui;
epis = jzero;
% find plasmonic root
trustr = realmax;
%kperp = sqrt(sphr.k.^2.*sphr.ep.*sphr.mu - sphr.beta.^2); %alpha c
kperp = sphr.k * sqrt(sphr.ep);
%proot = [];

%roots = apm(@sphrdispepinewtTM, 0, pi/2, 1, 0, 2, proot, 1e-5, sphr);
%roots = sort(roots);


% return two roots
%if(length(roots)>0)
if(sphr.orders ==0)
    alpha = pi/4*exp(-sphr.k*sphr.a/2);
else
     alpha = pi/2*exp(-sphr.k*sphr.a/2);
end
if(sphr.a*sphr.k<=1)
    proot(1:2) = newton(@sphrdispepinewtTM, [exp((pi-alpha)*1i)*sphr.k*sphr.a,exp(-alpha*1i)*sphr.k*sphr.a], trustr, sphr); %roots(1:2)
else
    proot(1:2) = newton(@sphrdispepinewtTM, [exp((pi-alpha)*1i)*3,exp(-alpha*1i)*3], trustr, sphr); %roots(1:2)
end
    %roots(2) = newton(@sphrdispepinewtTM, , trustr, sphr);
%end




% integration contours: avoid plasmonic root if necessary
[rad, cen] = contours(epis, proot);

% recalculate number of contours necessary
ncntrs = nroots; %ceil((nroots-length(proot)));
droots = zeros([1 ncntrs]);
droots(1) = proot(2);

% root search by contours centered on singularities
for k = 2:ncntrs
  droots(k) = search(sphr, epis,  droots(1:k-1), cen(k), rad(k));
end

% check for reality
roots = droots; %[proot droots];
if all(abs(imag(roots)./real(roots)) < eps)
  roots = real(roots);
end

% sort roots
[~, order] = sort(real(roots));
roots = roots(order);

% exclude unwanted roots
roots = roots(1:nroots);

function [rad, cen] = contours(sings, proot)
% lower boundary is 3/4 distance to next lower singularity 
dz = diff([0 abs(sings)]);
  

cen =  sings;
rad = dz/3;

function roots = search(sphr, sing, proot, cen, rad)
% course contour parameters
crstol = 1e-5;
conbuffer = 2;

% fine contour parameters
fintol = 1e-10;
sdistlim = 1e3*crstol;
fconsize = 1.5;

% newton-raphson parameters
trustr = 5;

 % detect if known roots are close to contour
% use only after first contour
if length(proot) > 1
  nearc = abs(abs(cen-proot)-rad) < conbuffer;

  % adjust contour
  if any(nearc)
    % existing contour
    llim = cen-rad; ulim = cen+rad;

    % new contour bisects root and next lowest contour
    buffroot = min(proot(nearc));
    lroot = max(sing(sing < buffroot));
    llim = (buffroot+lroot)/2;
    cen = (ulim+llim)/2; rad = (ulim-llim)/2;
  end
end

% detect singularities in contour
incs = abs(cen - sing) < rad;

sings = sing(incs);

%detect previous roots in large contour
incp = abs(cen - proot) < rad;
proots = proot(incp);

% large course contour to include at least two roots
phi = 0; obl = 1;
roots = apm(@sphrdispepinewtTM, cen, rad, obl, phi, sings, proots, crstol, sphr);
roots = sort(roots);


% return two roots
if(length(roots)>0)
roots(1) = newton(@sphrdispepinewtTM, roots(1), trustr, sphr); %roots(1:2);
else
    roots
end

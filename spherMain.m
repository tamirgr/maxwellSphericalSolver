function []=spherMain(firstL)
% define a sphrinder
sphr =SphereGeometry;

% background properties
sphr.ep = 1; %sphr.mu = 1.0;

% field properties
sphr.k = 1.0; %sphr.beta = 0.0;

% sphrinder properties
sphr.a = 1.0;
%sphr.mui = 1;

% sphrinder coordinates
sphr.x = 0.0; sphr.y = 0.0;

% orders of modes
morders = 0.0+firstL:10+firstL; 
nmodes = 20; % radial

% initialize mode orders
sphrl(1:length(morders)*nmodes) = sphr;
% order modes by azimuthal order in groups of radial order
mvar = num2cell(repmat(morders, [1 nmodes]));
[sphrl.orders] = deal(mvar{:});

% calculate sphrinder eigenvalues
epis = zeros([length(morders) nmodes]);
for k = 1:length(morders)
   if(firstL==0)
        epis(k,:) = disprootsepi3(sphrl(k), nmodes);
   else
        epis(k,:) = disprootsepi2(sphrl(k), nmodes);
   end
end


% % get modal fields
% for k = 1:length(sphrl)
%   sphrl(k).epi = epis(k);
%   sphrl(k) = modefield(sphrl(k));
% end
% 
% % calculate total fields
% npts = 100;
% lim = 2;
% [x, y, z] = meshgrid(linspace(-lim,lim,npts), linspace(-lim,lim,npts),linspace(-lim,lim,npts));
% 
% % calculate basis fields
% plotmode = 1;
% [E , H] = fieldpt2(x, y, z, sphrl(plotmode));
% 
% plot
pcolor(x, y, imag(Ezm))
viscircles([0 0], sphr.a, 'EdgeColor', 'k', 'LineWidth', 1, 'DrawBackgroundCircle', false);
shading flat

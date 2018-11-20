function [deltaR1,deltaT1,deltaPh1,deltaR2,deltaT2,deltaPh2] = checkBoundaryConditions(sphr,m,maxL)
    %direction r
    % define a sphrinder
if nargin < 3
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
    sphr.z = 0.0;
    m=3;
    maxL = 7;
end
% orders of modes
morders = 0.0:maxL; 
nmodes = 20; % radial

    % initialize mode orders
    sphrl(1:length(morders)*nmodes) = sphr;
% order modes by azimuthal order in groups of radial order
    mvar = num2cell(repmat(morders, [1 nmodes]));
    [sphrl.orders] = deal(mvar{:});

% calculate sphrinder eigenvalues
    epis1 = zeros([length(morders) nmodes]);
    epis2 = zeros([length(morders) nmodes]);
    for k = 1:length(morders)
          epis1(k,:) = disprootsepi2(sphrl(k), nmodes);
          epis2(k,:) = disprootsepi3(sphrl(k), nmodes);
    end

    deltaR1 = zeros(size(epis1));
    deltaT1 = zeros(size(epis1));
    deltaPh1 = zeros(size(epis1));
    deltaR2 = zeros(size(epis2));
    deltaT2 = zeros(size(epis2));
    deltaPh2 = zeros(size(epis2));
    outCoef = 1i/sphr.k;
    inCoef = 1i*sphr.ep./epis1;
    [r_in,theta_in,phi_in] = curlFXlmPoint(sphr.a,0.5,pi/4,epis1,sphrl(1));
    [r_out,theta_out,phi_out] = curlFXlmPointHenkel(sphr.a,0.5,pi/4,epis1,sphrl(1));
    deltaR1(:,:) = r_out(:,:,m).*outCoef - r_in(:,:,m).*inCoef;
    deltaT1(:,:) = theta_out(:,:,m).*outCoef-theta_in(:,:,m).*inCoef;
    deltaPh1(:,:) = phi_out(:,:,m).*outCoef-phi_in(:,:,m).*inCoef;
    
    [r_in,theta_in,phi_in] = curlFXlmPoint(sphr.a,0.5,pi/4,epis2,sphrl(1));
    [r_out,theta_out,phi_out] = curlFXlmPointHenkel(sphr.a,0.5,pi/4,epis2,sphrl(1));
    deltaR2(:,:) = r_out(:,:,m)*outCoef-r_in(:,:,m).*inCoef;
    deltaT2(:,:) = theta_out(:,:,m).*outCoef-theta_in(:,:,m).*inCoef;
    deltaPh2(:,:) = phi_out(:,:,m).*outCoef-phi_in(:,:,m).*inCoef;
end
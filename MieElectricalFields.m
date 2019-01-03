function [Ex,Ey,Ez,Hx,Hy,Hz] = MieElectricalFields(x, y, z,indx1,indx2,indx3,E0r,E0th,E0phi,H0r,H0th,H0phi)
if nargin < 3
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
        range =2;
        sphr.x = 0.0; sphr.y = 0.0;
        sphr.z = 0.0;
        x = linspace(-range,range,len);
        y = linspace(-range,range,len);
        z = linspace(-range,range,len);
        indx1 = 10; %l
        indx2 = 2; %n
        indx3 = 1; %m
        type = 1;
        nmodes = 20;
        %sphr.beta = 0.5;
        
        sphr.ordersN = nmodes;
        epiNL1 = zeros(indx1,nmodes);
        epiNL2 = zeros(indx1,nmodes);
         for k = 1:indx1
             sphr.orders = k;
             epiNL1(k,:) = disprootsepi2(sphr, nmodes);
             epiNL2(k,:) = disprootsepi3(sphr, nmodes);
         end
        sphr.orders = indx1; 
        mue = sqrt(1.5);
        
        [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
        
    E0r=ones(size(r));
        E0th=ones(size(th));
        E0phi=ones(size(phi));
    H0r=ones(size(r));
        H0th=ones(size(th));
        H0phi=ones(size(phi));
    else
    [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
end
    th = pi/2 -th;
   EiR = zeros(size(r));
   EiTh = zeros(size(th));
   EiPhi = zeros(size(phi));
   
   HiR = zeros(size(r));
   HiTh = zeros(size(th));
   HiPhi = zeros(size(phi));
   
   ELR = zeros(size(r));
   ELTh = zeros(size(th));
   ELPhi = zeros(size(phi));
   
   HLR = zeros(size(r));
   HLTh = zeros(size(th));
   HLPhi = zeros(size(phi));
   
   ESR = zeros(size(r));
   ESTh = zeros(size(th));
   ESPhi = zeros(size(phi));
   
   HSR = zeros(size(r));
   HSTh = zeros(size(th));
   HSPhi = zeros(size(phi));
   
   [a,b,c,d] = calcCoeff(sphr,epiNL1,indx2,indx1); 
   for l = 1:indx1
        [EiR,EiTh,EiPhi] = EiField(r,th,phi,sphr,epiNL1,indx2,l,indx3,EiR,EiTh,EiPhi,E0r,E0th,E0phi,ones(1,indx1),ones(1,indx1));
        [HiR,HiTh,HiPhi] = HiField(r,th,phi,sphr,epiNL1,indx2,l,indx3,HiR,HiTh,HiPhi,H0r,H0th,H0phi,ones(1,indx1),ones(1,indx1));
  
        [ELR,ELTh,ELPhi] = EiField(r,th,phi,sphr,epiNL1,indx2,l,indx3,ELR,ELTh,ELPhi,E0r,E0th,E0phi,c,d);
        [HLR,HLTh,HLPhi] = HiField(r,th,phi,sphr,epiNL1,indx2,l,indx3,HLR,HLTh,HLPhi,H0r,H0th,H0phi,d,c);
        
        [ESR,ESTh,ESPhi] = EiField(r,th,phi,sphr,epiNL1,indx2,l,indx3,ESR,ESTh,ESPhi,E0r,E0th,E0phi,-b,-a);
        [HSR,HSTh,HSPhi] = HiField(r,th,phi,sphr,epiNL1,indx2,l,indx3,HSR,HSTh,HSPhi,H0r,H0th,H0phi,-a,-b);
       
   end
   ext = r>sphr.a;
   [Ex,Ey,Ez] = mySph2cart((ESR + EiR).*(1-ext) + ext.*ELR,(ESTh + EiTh).*(1-ext) + ext.*ELTh,(ESPHi + EiPHi).*(1-ext) + ext.*ELPHi,th,phi);
   [Hx,Hy,Hz] = mySph2cart((HSR + HiR).*(1-ext) + ext.*HLR,(HSTh + HiTh).*(1-ext) + ext.*HLTh,(HSPHi + HiPHi).*(1-ext) + ext.*HLPHi,th,phi);
   
    ExR = real(Ex);
     EyR = imag(Ey);
     EzR = real(Ez);
     
    figure(1);
    slice(X,Y,Z,ExR,x(len/2),x(len/2),x(len/2)+0.1);
    %slice(X,Y,Z,ExR,dispx,dispx,dispx);
    colorbar();
    xlabel('x');
    ylabel('y');
    zlabel('z');
    %caxis([-0.1,0.1]);
    shading interp
    title(sprintf('ExReal n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(1,sprintf('sliceExRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(1,sprintf('sliceExRealn%dl%dm%d.fig',indx2,indx1,indx3));
    
     
    figure(2);
    slice(X,Y,Z,EyR,x(len/2),x(len/2),x(len/2));
    %slice(X,Y,Z,EyR,dispx,dispx,x(len/2));
    shading interp
    colorbar();
        xlabel('x');
    ylabel('y');
    zlabel('z');
    caxis([-0.1,0.1]);
    title(sprintf('EyReal n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(2,sprintf('sliceEyRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(2,sprintf('sliceEyRealn%dl%dm%d.fig',indx2,indx1,indx3));
    %slice(X,Y,Z,Ey,[-1,0,1],0,-1);
    figure(3);
    slice(X,Y,Z,EzR,x(len/2),x(len/2),x(len/2));
    %slice(X,Y,Z,EzR,dispx,dispx,x(len/2));
        xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar();
    caxis([-0.1,0.1]);
    shading interp
    title(sprintf('EzReal n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(3,sprintf('sliceEzRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(3,sprintf('sliceEzRealn%dl%dm%d.fig',indx2,indx1,indx3));

end
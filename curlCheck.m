function [res] = curlCheck(x, y, z, epiNL1, sphr,indx1,indx2,indx3,type)
    if nargin < 3
        sphr =SphereGeometry;

        % background properties
        sphr.ep = 1; %sphr.mu = 1.0;

        % field properties
        sphr.k = 0.1; %sphr.beta = 0.0;

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
        indx1 = 1; %l
        indx2 = 2; %n
        indx3 = 1; %m
        type = 1;
        nmodes = 20;
        %sphr.beta = 0.5;
        
        sphr.ordersN = nmodes;
        epiNL1 = zeros(indx1,nmodes);
    %    epiNL2 = zeros(indx1,nmodes);
         for k = 1:10
             sphr.orders = k-1;
             epiNL1(k,:) = disprootsepi2(sphr, nmodes);
    %         epiNL2(k,:) = disprootsepi3(sphr, nmodes);
         end
        sphr.orders = indx1; 
    end
    [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
   th = pi/2 -th;

    res = zeros(length(r),length(th),length(phi),3);

    
    dispx = [x(floor(len/3)),x(floor(len/3*2))];
    f = Fnl(r,sphr,epiNL1,indx2,indx1);

        resR = zeros(size(r));
        [xlmTh,xlmPhi] = Xlm(th,phi,indx1,indx3);
        

            resTh = xlmTh.*f;
            resPhi = xlmPhi.*f;

     [R,Th,Phi] = TEField(r,th,phi,sphr,epiNL1,indx2,indx1,indx3);
    [Ex, Ey, Ez] =  mySph2cart(R,Th,Phi,th,phi);

    %[Ex,Ey,Ez] = mySph2cart(resR,resTh,resPhi,th,phi);
    
     ExR = real(Ex);
     EyR = real(Ey);
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
    title(sprintf('ExReal n=%d, l=%d, m=%d\n eigen mode  = %0.2f + %0.2f*i ',indx2,indx1,indx3,real(epiNL1(indx1+1,indx2)),imag(epiNL1(indx1+1,indx2))));
    saveas(1,sprintf('sliceExRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(1,sprintf('sliceExRealn%dl%dm%d.fig',indx2,indx1,indx3));
    
     
    figure(2);
    slice(X,Y,Z,EyR,x(len/2),x(len/2),x(len/2)+0.05);
    %slice(X,Y,Z,EyR,dispx,dispx,x(len/2));
    shading interp
    colorbar();
        xlabel('x');
    ylabel('y');
    zlabel('z');
 %   caxis([-0.1,0.1]);
    title(sprintf('EyReal n=%d, l=%d, m=%d\n eigen mode = %0.2f + %0.2f*i ',indx2,indx1,indx3,real(epiNL1(indx1+1,indx2)),imag(epiNL1(indx1+1,indx2))));
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
  %  caxis([-0.1,0.1]);
    shading interp
    title(sprintf('EzReal n=%d, l=%d, m=%d\n eigen mode = %0.2f + %0.2f*i ',indx2,indx1,indx3,real(epiNL1(indx1+1,indx2)),imag(epiNL1(indx1+1,indx2))));
    saveas(3,sprintf('sliceEzRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(3,sprintf('sliceEzRealn%dl%dm%d.fig',indx2,indx1,indx3));
     
%     ext = (r<sphr.a);
%     [resRIn,resThIn,resPhiIn] = curlFXlmPoint1(r,th,phi,epiNL2,sphr,indx2,indx1,indx3);
%     [resROut,resThOut,resPhiOut] = curlFXlmPointHenkel2(r,th,phi,epiNL2,sphr,indx2,indx1,indx3);    
%     resROut = 1i/sphr.k*(ext.*resRIn/epiNL2(indx1+1,indx2)*sphr.ep + (~ext).*resROut);
%     resThOut = 1i/sphr.k*(ext.*resThIn/epiNL2(indx1+1,indx2)*sphr.ep + (~ext).*resThOut);
%     resPhiOut = 1i/sphr.k*(ext.*resPhiIn/epiNL2(indx1+1,indx2)*sphr.ep + (~ext).*resPhiOut);
%     %[Ex,Ey,Ez] = mySph2cart(resROut,resThOut,resPhiOut,th,phi);
%     [ExRot,EyRot,EzRot] = mySph2cart(resROut,resThOut,resPhiOut,th,phi);
%     
%     [ExR,EyR,EzR] = mySph2cart(resR,resTh,resPhi,th,phi);
%     
%      ExR = real(ExRot);
%      EyR = real(EyRot);
%      EzR = real(EzRot);
%      
%     % EyR = real(resPhi);
%     % EzR = real(resTh);
%     %ExI = imag(Ex);
%     %EyI = imag(Ey);
%     %EzI = imag(Ez);
%     %quiver3(x,y,z,x,y,z);
%     %[X,Y,Z] = meshgrid(x,y,z);
%     
%     
%     figure(4);
%     slice(X,Y,Z,ExR,dispx,dispx,dispx);
%     colorbar();
%     shading interp
%     title(sprintf('ExReal n=%d, l=%d, m=%d\n eigen mode = %0.2f + %0.2f*i ',indx2,indx1,indx3,real(epiNL2(indx1+1,indx2)),imag(epiNL2(indx1+1,indx2))));
%     saveas(4,sprintf('sliceBxRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
%     saveas(4,sprintf('sliceBxRealn%dl%dm%d.fig',indx2,indx1,indx3));
%     
%      
%     figure(5);
%     slice(X,Y,Z,EyR,dispx,dispx,x(len/2));
%     shading interp
%     colorbar();
%     title(sprintf('EyReal n=%d, l=%d, m=%d\n eigen mode = %0.2f + %0.2f*i ',indx2,indx1,indx3,real(epiNL2(indx1+1,indx2)),imag(epiNL2(indx1+1,indx2))));
%     saveas(5,sprintf('sliceByRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
%     saveas(5,sprintf('sliceByRealn%dl%dm%d.fig',indx2,indx1,indx3));
%     %slice(X,Y,Z,Ey,[-1,0,1],0,-1);
%     figure(6);
%     slice(X,Y,Z,EzR,dispx,dispx,x(len/2));
%     colorbar();
%     shading interp
%     title(sprintf('EzReal n=%d, l=%d, m=%d\n eigen mode = %0.2f + %0.2f*i ',indx2,indx1,indx3,real(epiNL2(indx1+1,indx2)),imag(epiNL2(indx1+1,indx2))));
%     saveas(6,sprintf('sliceBzRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
%     saveas(6,sprintf('sliceBzRealn%dl%dm%d.fig',indx2,indx1,indx3));
%     %slice(X,Y,Z,Ez,[-1,0,1],0,-1);
    

end
function [res] = plotField(x, y, z, epiNL1, sphr,indx1,indx2,indx3,type)
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
        [x,y,z]= sphere(len);
        indx1 = 2; %l
        indx2 = 2; %n
        indx3 = 1; %m
        type = 1;
        nmodes = 20;
        %sphr.beta = 0.5;
        
        sphr.ordersN = nmodes;
        epiNL1 = zeros(indx1,nmodes);
        epiNL2 = zeros(indx1,nmodes);
         for k = 1:10
             sphr.orders = k;
             epiNL1(k,:) = disprootsepi2(sphr, nmodes);
             epiNL2(k,:) = disprootsepi3(sphr, nmodes);
         end
        sphr.orders = indx1; 
    end
    %[X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(x,y,z);
   th = pi/2 -th;
    if(type == 1)
        res = zeros(length(r),length(th),length(phi),3);
    else
         res = zeros([size(epiNL1),2*size(epiNL1,1)]);
    end
    
    %dispx = [x(floor(len/3)),x(floor(len/3*2))];
    f = Fnl(r,sphr,epiNL1,indx2,indx1);
    if(type == 1)
        resR = zeros(size(r));
        [xlmTh,xlmPhi] = Xlm(th,phi,indx1,indx3);
        
       % for r1 = 1:length(r)
            resTh = xlmTh.*f;
            resPhi = xlmPhi.*f;
       % end
    else
        for l = 1:size(epiNL1,1)
            value = Ylm(indx1,th,phi);
            res(l,:,:) = bsxfun(@times,f(l,:,r(indx1)),value(indx2,indx3,:));
        end
    end
    %[Ex,Ey,Ez] = mySph2cart(resR,resTh,resPhi,th,phi);
    
     %ExR = real(Ex);
     %EyR = real(Ey);
     %EzR = real(Ez);
     
    figure(1);
    surf(x,y,z,real(resTh));
    %slice(X,Y,Z,ExR,x(len/2),x(len/2),x(len/2));
    %slice(X,Y,Z,ExR,dispx,dispx,dispx);
    colorbar();

    caxis([-0.1,0.1]);
    shading interp
    title(sprintf('sphereTh n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(1,sprintf('sphereThRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(1,sprintf('sphereThRealn%dl%dm%d.fig',indx2,indx1,indx3));
    
     
    figure(2);
    surf(x,y,z,real(resPhi));
    %slice(X,Y,Z,EyR,x(len/2),x(len/2),x(len/2));
    %slice(X,Y,Z,EyR,dispx,dispx,x(len/2));
    shading interp
    colorbar();
    caxis([-0.1,0.1]);
    title(sprintf('Phi n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(2,sprintf('spherePhiRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(2,sprintf('spherePhiRealn%dl%dm%d.fig',indx2,indx1,indx3));
    %slice(X,Y,Z,Ey,[-1,0,1],0,-1);
     
    ext = (r<1);
    [resRIn,resThIn,resPhiIn] = curlFXlmPoint1(r,th,phi,epiNL2,sphr,indx2,indx1,indx3);
    [resROut,resThOut,resPhiOut] = curlFXlmPointHenkel1(r,th,phi,epiNL2,sphr,indx2,indx1,indx3);    
    resROut = 1i/sphr.k*(ext.*resRIn/epiNL2(indx1,indx2)*sphr.ep + (~ext).*resROut);
    resThOut = 1i/sphr.k*(ext.*resThIn/epiNL2(indx1,indx2)*sphr.ep + (~ext).*resThOut);
    resPhiOut = 1i/sphr.k*(ext.*resPhiIn/epiNL2(indx1,indx2)*sphr.ep + (~ext).*resPhiOut);
    %[Ex,Ey,Ez] = mySph2cart(resROut,resThOut,resPhiOut,th,phi);
    %[ExRot,EyRot,EzRot] = mySph2cart(resROut,resThOut,resPhiOut,th,phi);
    
    
     figure(3);
    surf(x,y,z,real(resThOut));
    %slice(X,Y,Z,ExR,x(len/2),x(len/2),x(len/2));
    %slice(X,Y,Z,ExR,dispx,dispx,dispx);
    colorbar();

    caxis([-0.1,0.1]);
    shading interp
    title(sprintf('curlsphereTh n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(1,sprintf('sphereThRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(1,sprintf('sphereThRealn%dl%dm%d.fig',indx2,indx1,indx3));
    
     
    figure(4);
    surf(x,y,z,real(resPhiOut));
    %slice(X,Y,Z,EyR,x(len/2),x(len/2),x(len/2));
    %slice(X,Y,Z,EyR,dispx,dispx,x(len/2));
    shading interp
    colorbar();
    caxis([-0.1,0.1]);
    title(sprintf('curlPhi n=%d, l=%d, m=%d',indx2,indx1,indx3));
    saveas(2,sprintf('curlspherePhiRealn%dl%dm%d.jpeg',indx2,indx1,indx3));
    saveas(2,sprintf('curlspherePhiRealn%dl%dm%d.fig',indx2,indx1,indx3));

end
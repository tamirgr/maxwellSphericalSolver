function  printCurlQuivers(x, y, z, epiNL1, sphr,indx1,indx2,indx3,type)
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
        indx1 = 1; %l
        indx2 = 1; %n
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
    [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
   th = pi/2 -th;
    if(type == 1)
        res = zeros(length(r),length(th),length(phi),3);
    else
         res = zeros([size(epiNL1),2*size(epiNL1,1)]);
    end
    
    dispx = [x(floor(len/3)),x(floor(len/3*2))];
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
    [Ex,Ey,Ez] = mySph2cart(resR,resTh,resPhi,th,phi);
    
     ExR = real(Ex);
     EyR = real(Ey);
     EzR = real(Ez);
     
    ext = (r<sphr.a);
    ext1 = r<sphr.a*1.05-ext;
    ext2 = ext>(r<sphr.a*0.95);
    
    [resRIn,resThIn,resPhiIn] = curlFXlmPoint2(r,th,phi,epiNL2,sphr,indx2,indx1,indx3);
    [resRO,resThO,resPhiO] = curlFXlmPointHenkel2(r,th,phi,epiNL2,sphr,indx2,indx1,indx3);    
    resROut = 1i/sphr.k*(ext2.*resRIn/epiNL2(indx1,indx2)*sphr.ep );
    resThOut = 1i/sphr.k*(ext2.*resThIn/epiNL2(indx1,indx2)*sphr.ep);
    resPhiOut = 1i/sphr.k*(ext2.*resPhiIn/epiNL2(indx1,indx2)*sphr.ep);
    [ExRot,EyRot,EzRot] = mySph2cart(resROut.*ext2,resThOut,resPhiOut,th,phi);
    ExR = real(ExRot);
    EyR = real(EyRot);
    EzR = real(EzRot);

    quiver3(X,Y,Z,ExR,EyR,EzR);
    title('in');
    resROut = 1i/sphr.k*(ext1.*resRO);
    resThOut = 1i/sphr.k*(ext1.*resThO);
    resPhiOut = 1i/sphr.k*(ext1.*resPhiO);
    [ExRot,EyRot,EzRot] = mySph2cart(resROut.*ext2,resThOut,resPhiOut,th,phi);


     ExR = real(ExRot);
     EyR = real(EyRot);
     EzR = real(EzRot);
    figure;
    quiver3(X,Y,Z,ExR,EyR,EzR);
    title('out');


end
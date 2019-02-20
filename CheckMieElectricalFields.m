function [diffETh,diffEphi,diffHTh,diffHphi] = CheckMieElectricalFields(x, y, z,indx1,indx2,indx3,E0r,E0th,E0phi,H0r,H0th,H0phi)

if nargin < 3
        sphr =SphereGeometry;
        
        % background properties
        sphr.ep = 1; %sphr.mu = 1.0;

        % field properties
        sphr.k = 5.0; %sphr.beta = 0.0;

        % sphrinder properties
        sphr.a = 1.0;
        %sphr.mui = 1;
        %sphr.mu = 2;

        % sphrinder coordinates
        len = 40;
        range =2;
        sphr.x = 0.0; sphr.y = 0.0;
        sphr.z = 0.0;
        x = linspace(-range,range,len);
        y = linspace(-range,range,len);
        z = linspace(-range,range,len);
        indx1 = 2; %l
        indx2 = 10; %n
        indx3 = 1; %m
        type = 1;
        nmodes = 20;
        %sphr.beta = 0.5;
        
        sphr.ordersN = nmodes;
%         epiNL1 = zeros(indx1,nmodes);
%         epiNL2 = zeros(indx1,nmodes);
%          for k = 1:indx1+1
%              sphr.orders = k;
%              epiNL1(k,:) = disprootsepi2(sphr, nmodes);
%              epiNL2(k,:) = disprootsepi3(sphr, nmodes);
%          end
        sphr.orders = indx1; 
        mu = ones(indx1+1,indx2)*(sqrt(1.5));
        
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
 mu = ones(indx1+1,indx2)*(sqrt(1.5));
    close all;
th = pi/2 -th;
   EiR = zeros(size(r));
   EiTh = zeros(size(th));
   EiPhi = zeros(size(phi));
   r = sphr.a*ones(size(r)); 
   %HiR = zeros(size(r));
  % HiTh = zeros(size(th));
  % HiPhi = zeros(size(phi));
   
   ELR = zeros(size(r));
   ELTh = zeros(size(th));
   ELPhi = zeros(size(phi));
   
   HLR = zeros(size(r));
   HLTh = zeros(size(th));
   HLPhi = zeros(size(phi));
   
   ESR = zeros(size(r));
   ESTh = zeros(size(th));
   ESPhi = zeros(size(phi));
   
   %HSR = zeros(size(r));
   %HSTh = zeros(size(th));
   %HSPhi = zeros(size(phi));
   
   [a,b,c,d] = calcCoeff(sphr,mu,indx1); 
   %for l = 1:indx1
    l = indx1;
    n = indx2;
   % for n = 1:indx2
%         EiR = exp(1i*sphr.k*r.*cos(th)).*sin(th).*cos(phi);
%         EiTh = exp(1i*sphr.k*r.*cos(th)).*cos(th).*cos(phi);
%         EiPhi = -exp(1i*sphr.k*r.*cos(th)).*sin(phi);
%         
     %   HiR = 1i*sphr.k *exp(1i*sphr.k*r.*sin(th)).*sin(th).*sin(phi);
     %   HiTh = -1i*sphr.k*exp(1i*sphr.k*r.*sin(th)).*cos(th).*sin(phi);
     %   HiPhi = 1i*sphr.k *exp(1i*sphr.k*r.*sin(th)).*cos(phi);
   
        [Mr,Mth,Mphi] = MOField(r,th,phi,sphr,mu(1,1),l,indx3);
        [Nr,Nth,Nphi] = NEField(r,th,phi,sphr,mu(1,1),l,indx3);
        [EiR,EiTh,EiPhi] = EiField(l,EiR,EiTh,EiPhi,E0r,E0th,E0phi,ones(1,indx1),ones(1,indx1),Mr,Mth,Mphi,Nr,Nth,Nphi);
        
        [ELR,ELTh,ELPhi] = EiField(l,ELR,ELTh,ELPhi,E0r,E0th,E0phi,c,d,Mr,Mth,Mphi,Nr,Nth,Nphi);
        
        [ESR,ESTh,ESPhi] = EiField(l,ESR,ESTh,ESPhi,E0r,E0th,E0phi,-b,-a,Mr,Mth,Mphi,Nr,Nth,Nphi);
    
%         [Mr,Mth,Mphi] = MEField(r,th,phi,sphr,l,indx3);
%         [Nr,Nth,Nphi] = NOField(r,th,phi,sphr,l,indx3);
        
          % [HiR,HiTh,HiPhi] = HiField(sphr,epiNL1,indx2,l,HiR,HiTh,HiPhi,H0r,H0th,H0phi,ones(1,indx1),ones(1,indx1),Mr,Mth,Mphi,Nr,Nth,Nphi);
      %[HLR,HLTh,HLPhi] = HiField(sphr,mu,indx2,l,HLR,HLTh,HLPhi,H0r,H0th,H0phi,d,c,Mr,Mth,Mphi,Nr,Nth,Nphi);
    
       %   [HSR,HSTh,HSPhi] = HiField(sphr,mu,indx2,l,HSR,HSTh,HSPhi,H0r,H0th,H0phi,-a,-b,Mr,Mth,Mphi,Nr,Nth,Nphi);
       
   %end
 %  ext = r>sphr.a;
  %[Ex,Ey,Ez] = mySph2cart( EiR ,(EiTh ) ,(EiPhi ) ,th,phi);
  
  % [Ex,Ey,Ez] = mySph2cart((ESR ).*(~ext) + ext.*ELR,(ESTh ).*(~ext) + ext.*ELTh,(ESPhi ).*(~ext) + ext.*ELPhi,th,phi);
%   [Hx,Hy,Hz] = mySph2cart((HSR + HiR).*(1-ext) + ext.*HLR,(HSTh + HiTh).*(1-ext) + ext.*HLTh,(HSPhi + HiPhi).*(1-ext) + ext.*HLPhi,th,phi);
   diffETh = max(max(max(abs((EiTh + ESTh - ELTh)./ELTh)))) ;
   diffEphi = max(max(max(abs((EiPhi + ESPhi - ELPhi)./ELPhi)) ));
   
   epsback = 1.0;
    close all
    L=l;
    N=n;
    epiNL = zeros(L,N);
    for l=1:L+1
        sphr.orders = l-1;
        epiNL(l,:) = disprootsepi2(sphr, N);
    end

  %  [ExRot,EyRot,EzRot] = genTEField2(epsback, sphr, n, l, m, len, range);
    [ERIn, EThIn, EPhiIn] = TEFieldDebug(r,th,phi,sphr,epiNL,n,L,1,0);
    [ERout, EThout, EPhiout] = TEFieldDebug(r,th,phi,sphr,epiNL,n,L,1,1);
                
   
   diffHTh = max(max(max( (EThIn-EThout)./EThIn)));
   diffHphi = max(max(max((EPhiIn-EPhiout)./EPhiIn)));

end
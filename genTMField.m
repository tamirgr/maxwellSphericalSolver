function [ Ex, Ey, Ez ] = genTMField( ExZero, EyZero, EzZero, Jx, Jy, Jz, epsinc, epsback, sphr, N, L)
%genTEField - calculates the total field as described in the paper 'Generalizing normal mode expansion of electromagnetic Greens tensor to open systems'
%   Ex, Ey, Ez - the eigenmodes of the system. |Em>
%   eigenvals - the corresponding eigenvalues of the system.
%   J - the source of the field. |J>
%   epsinc - the inclusion permittivity of the material
%   epsback - the background permittivity 
    %% Init
    len = 60;
    range = 2;

    if nargin < 6
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
        sphr.x = 0.0; sphr.y = 0.0;
        sphr.z = 0.0;

        L = 10;
        N = 2;

        sphr.ordersN = N;
        
        ExZero = zeros(len,len,len);
        EyZero = zeros(len,len,len);
        EzZero = zeros(len,len,len);
        Jx = ones(len,len,len);
        Jy = ones(len,len,len);
        Jz = ones(len,len,len);
        epsinc = sqrt(1.5);
        epsback = 1.0;
    end
    
    x = linspace(-range,range,len);
    y = linspace(-range,range,len);
    z = linspace(-range,range,len);
    [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
    th = pi/2 - th;
    
    %% Debug Parameters (Delete section when finished
     m=2;
    
    %% Calculate Eigenmodes and Eigenvalues

    epiNL = zeros(L,sphr.ordersN);
    for l=1:L
        sphr.orders = l;
        epiNL(l,:) = disprootsepi3(sphr, sphr.ordersN);
    end
    
    EmTh = X*0;
    EmPhi = Y*0;
    EmR = Z*0;
    
    n = N;
%     for n=1:N
        for l=1:L
%             for m=-l:l
                [ER,ETh, EPhi] = TMField(r,th,phi,sphr,epiNL,n,l,m);
                EmTh  = EmTh  + (epsinc-epsback)/(epiNL(l,n)-epsback)/(epiNL(l,n)-epsinc)*ETh;
                EmPhi = EmPhi + (epsinc-epsback)/(epiNL(l,n)-epsback)/(epiNL(l,n)-epsinc)*EPhi;
                EmR   = EmR   + (epsinc-epsback)/(epiNL(l,n)-epsback)/(epiNL(l,n)-epsinc)*ER;
%             end
        end
%     end
    
    %% Summation over all eigenmodes
   
	[Emx, Emy, Emz] =  mySph2cart(EmR,EmTh,EmPhi,th,phi);

	emBraJx = permute(conj(Emx),[2,1,3]).*Jx; % <Em|J> x-direction
    emBraJy = permute(conj(Emy),[2,1,3]).*Jy; % <Em|J> x-direction
    emBraJz = permute(conj(Emz),[2,1,3]).*Jz; % <Em|J> x-direction
    Ex = ExZero + 1i./sphr.k.*Emx.*emBraJx; %|Em><Em|J> x-direction
    Ey = EyZero + 1i./sphr.k.*Emy.*emBraJy; %|Em><Em|J> y-direction
    Ez = EzZero + 1i./sphr.k.*Emz.*emBraJz; %|Em><Em|J> z-direction
end


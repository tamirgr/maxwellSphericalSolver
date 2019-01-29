function [ Ex, Ey, Ez ] = genTMField2( epsback, sphr, N, L, M, len, range)
%genTEField - calculates the total field as described in the paper 'Generalizing normal mode expansion of electromagnetic Greens tensor to open systems'
% specifially - equation 21.
%   Ex, Ey, Ez - the eigenmodes of the system. |Em>
%   eigenvals - the corresponding eigenvalues of the system.
%   J - the source of the field. |J>
%   epsinc - the inclusion permittivity of the material
%   epsback - the background permittivity 
    %% Init

    if nargin < 7
        len = 40;
        range = 2;
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
        
        epsback = 1.0;
    else
    %% Debug Parameters (Delete section when finished
        m = M;
    end
    
    x = linspace(-range,range,len);
    y = linspace(-range,range,len);
    z = linspace(-range,range,len);
    [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
    th = pi/2 - th;
    internal = (r<=sphr.a);
    tmp = (ones(len,len,len)-internal)*sphr.a^3*4*pi/3;
    integmat = tmp + internal.*r.^3*4*pi/3;
    %% Calculate Eigenmodes and Eigenvalues

    epiNL = zeros(L,N);
    for l=1:L
        sphr.orders = l-1;
        epiNL(l,:) = disprootsepi3(sphr, N);
    end
    
    EmTh = X*0;
    EmPhi = Y*0;
    EmR = Z*0;
    l = 1; %debug line
    epsi0 = real(epiNL(l+1,1));
    [ErZero,EthZero, EphiZero] = TMField(r,th,phi,sphr,epsi0,1,l,1); %calculate an eigenmode seperately for E0
    [Era0,Etha0, Ephia0] = TMField(sphr.a,pi/2,0,sphr,epsi0,1,l,1); %calculate the value of the field on sphere
%     ErZero = real(ErZero);
%     EthZero = real(EthZero);
%     EphiZero = real(EphiZero);

    for n=1:N
%         for l=1:L
%             for m=-l:l
                [ER, ETh, EPhi] = TMField(r,th,phi,sphr,epiNL,n,l,m);
                epco = (epsi0-epsback)/(epiNL(l+1,n)-epsi0); 
                EtR = ER.*internal + (1-internal)*Era0;
                EtTh = ETh.*internal + (1-internal)*Etha0;
                EtPhi = EPhi.*internal + (1-internal)*Ephia0;
                ER = conj(EtR).*ErZero.*integmat;
                ETh = conj(EtTh).*EthZero.*integmat;
                EPhi = conj(EtPhi).*EphiZero.*integmat;
                EmR   = EmR   + epco*ER;
                EmTh  = EmTh  + epco*ETh;
                EmPhi = EmPhi + epco*EPhi;
%             end
%         end
    end
    
    %% Summation over all eigenmodes
   
	[Emx, Emy, Emz] =  mySph2cart(EmR,EmTh,EmPhi,th,phi);
    
    [ExZero,EyZero,EzZero] = mySph2cart(ErZero,EthZero,EphiZero,th,phi);
%     disp = [x(floor(len/3)),x(floor(len/3*2))];
    disp = [x(floor(len/2))];
    displayFields( real(ExZero) , real(EyZero) , real(EzZero) ,X,Y,Z, 1,1,1,disp);

    sizefac = 1;%0*10^4;
    sizefac0 = 1;%0*10^4;
    Ex = ExZero*sizefac0 + 1i./sphr.k.*Emx.*sizefac; %|Em><Em|J> x-direction
    Ey = EyZero*sizefac0 + 1i./sphr.k.*Emy.*sizefac; %|Em><Em|J> y-direction
    Ez = EzZero*sizefac0 + 1i./sphr.k.*Emz.*sizefac; %|Em><Em|J> z-direction
end
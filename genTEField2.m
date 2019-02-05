function [ Ex, Ey, Ez ] = genTEField2( epsback, sphr, N, L, M, len, range)
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
    %% Calculate Eigenmodes and Eigenvalues

    epiNL = zeros(L,N);
    for l=1:L+1
        sphr.orders = l-1;
        epiNL(l,:) = disprootsepi2(sphr, N);
    end
    
    Emx = X*0;
    Emy = Y*0;
    Emz = Z*0;
    l = L; %debug line
%     epsi0 = real(epiNL(l+1,1));
    epsi0 = epiNL(l+1,l)*0.9;
%     epsi0 = 20;
    [ErZero,EthZero, EphiZero] = TEField(r,th,phi,sphr,epsi0,1,l,1); %calculate an eigenmode seperately for E0
    [ExZero,EyZero,EzZero] = mySph2cart(ErZero,EthZero,EphiZero,th,phi);
    
    for n=1:N
%         for l=1:L
%             for m=-l:l
                [ER, ETh, EPhi] = TEField(r,th,phi,sphr,epiNL,n,l,m);
                [EX, EY, EZ] = mySph2cart(ER, ETh, EPhi,th,phi);
                
                epco = (epsi0-epsback)/(epiNL(l+1,n)-epsi0); % epsilon coefficient
                
                eigenco = integrateJ2R2(l, sphr.a, epsi0, epiNL(l+1,n));
                
                Emx = EX.*epco.*eigenco;
                Emy = EY.*epco.*eigenco;
                Emz = EZ.*epco.*eigenco;
%             end
%         end
    end
    
    %% Summation over all eigenmodes
   
% 	[Emx, Emy, Emz] =  mySph2cart(EmR,EmTh,EmPhi,th,phi);
    
%     [ExZero1,EyZero1,EzZero1] = mySph2cart(ErZero1,EthZero1,EphiZero1,th,phi);
%     disp = [x(floor(len/3)),x(floor(len/3*2))];
%     displayFields( real(ExZero1) , real(EyZero1) , real(EzZero1) ,X,Y,Z, 2,2,2,disp);
    [E1R,E1TH, E1PHI] = TEField(r,th,phi,sphr,epiNL(2,1),1,1,1); %calculate an eigenmode seperately for E0
    [Ex1,Ey1,Ez1] = mySph2cart(E1R,E1TH, E1PHI,th,phi);
    disp = [x(floor(len/2))];
    displayFields( real(Ex1) , real(Ey1) , real(Ez1) ,X,Y,Z, 1,1,1,disp,1);

    sizefac = 1;%0*10^4;
    sizefac0 = 1;%0*10^4;
    Ex = ExZero*sizefac0 + 1i./sphr.k.*Emx.*sizefac; %|Em><Em|J> x-direction
    Ey = EyZero*sizefac0 + 1i./sphr.k.*Emy.*sizefac; %|Em><Em|J> y-direction
    Ez = EzZero*sizefac0 + 1i./sphr.k.*Emz.*sizefac; %|Em><Em|J> z-direction
end
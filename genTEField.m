function [ Ex, Ey, Ez ] = genTEField( ExZero, EyZero, EzZero, Jx, Jy, Jz, epsinc, epsback, sphr, N, L)
%genTEField - calculates the total field as described in the paper 'Generalizing normal mode expansion of electromagnetic Greens tensor to open systems'
%   Ex, Ey, Ez - the eigenmodes of the system. |Em>
%   eigenvals - the corresponding eigenvalues of the system.
%   J - the source of the field. |J>
%   epsinc - the inclusion permittivity of the material
%   epsback - the background permittivity 

    %% Init
    len = 60;
    range = 2;
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
        epiNL(l,:) = disprootsepi2(sphr, sphr.ordersN);
    end
    
    EmTh = X*0;
    EmPhi = Y*0;
    
    for n=1:N
        for l=1:L
%             for m=-l:l
                [ETh, EPhi] = TEField(r,th,phi,sphr,epiNL,n,l,m);
                EmTh  = EmTh  + (epsinc-epsback)/(epiNL(l,n)-epsback)/(epiNL(l,n)-epsinc)*ETh;
                EmPhi = EmPhi + (epsinc-epsback)/(epiNL(l,n)-epsback)/(epiNL(l,n)-epsinc)*EPhi;
%             end
        end
    end
    
    %% Summation over all eigenmodes
   
	[Emx, Emy, Emz] =  mySph2cart(0,EmTh,EmPhi,th,phi);

	emBraJx = ctranspose(Emx)*Jx; % <Em|J> x-direction
    emBraJy = ctranspose(Emy)*Jy; % <Em|J> y-direction
    emBraJz = ctranspose(Emz)*Jz; % <Em|J> z-direction
    Ex = ExZero + 1i./sphr.k.*Emx.*emBraJx; %|Em><Em|J> x-direction
    Ey = EyZero + 1i./sphr.k.*Emy.*emBraJy; %|Em><Em|J> y-direction
    Ez = EzZero + 1i./sphr.k.*Emz.*emBraJz; %|Em><Em|J> z-direction
end


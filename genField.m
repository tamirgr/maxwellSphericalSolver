function [ Ex, Ey, Ez ] = genField( Emx, Emy, Emz, eigenvals, Jx, Jy, Jz, epsinc, epsback, k)
%genField - calculates the total field as described in the paper 'Generalizing normal mode expansion of electromagnetic Greens tensor to open systems'
%   eigenmodes - the eigenmodes of the system. |Em>
%   eigenvals - the corresponding eigenvalues of the system.
%   J - the source of the field. |J>
%   epsinc - the inclusion permittivity of the material
%   epsback - the background permittivity 

    lenv = length(eigenvals);
    lenx = length(Emx);
    leny = length(Emx);
    lenz = length(Emx);
    if lenv~=lenx || lenx~=leny || leny~=lenz 
        Ex = 0;
        Ey = 0;
        Ez = 0;
        return
    end
    emBraJx = ctranspose(Emx)*Jx; % <Em|J> x-direction
    emBraJy = ctranspose(Emy)*Jy; % <Em|J> y-direction
    emBraJz = ctranspose(Emz)*Jz; % <Em|J> z-direction
    eps = (epsinc-epsback)./(eigenmodes-epsback)./(eigenmodes-epsinc);
    Ex = 1i./k.*sum(eps.*Emx)*emBraJx; %|Em><Em|J> x-direction
    Ey = 1i./k.*sum(eps.*Emy)*emBraJy; %|Em><Em|J> y-direction
    Ez = 1i./k.*sum(eps.*Emz)*emBraJz; %|Em><Em|J> z-direction
end


function [ res ] = genField( eigenmodes, eigenvals, J, epsinc, epsback, k)
%genField - calculates the total field as described in the paper 'Generalizing normal mode expansion of electromagnetic Greens tensor to open systems'
%   eigenmodes - the eigenmodes of the system. |Em>
%   eigenvals - the corresponding eigenvalues of the system.
%   J - the source of the field. |J>
%   epsinc - the inclusion permittivity of the material
%   epsback - the background permittivity 

    lenv = length(eigenvals);
    lenm = length(eigenmodes);
    if lenv~=lenm
        res = 0;
        return
    end
    eigenmodesBra = ctranspose(eigenmodes); % <Em|
    eps = (epsinc-epsback)./(eigenmodes-epsback)./(eigenmodes-epsinc);
    res = 1i/k*sum(eps.*eigenmodes).*eigenmodesBra.*J; %|Em><Em|J>
end


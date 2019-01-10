function [ Jx, Jy, Jz] = genWave(len, range, wavetype, E0, Rs, k)
%genWave generates the matrices representing the desired wave
%   len - size of result matrices
%   range - real boundaries of matrices
%   wavetype - type of source wave
%   E0 - strength of wave
%   Rs - location of source

    %% Init
    Jx = ones(len,len,len);
    Jy = ones(len,len,len);
    Jz = ones(len,len,len);

    x = linspace(-range,range,len);
    y = linspace(-range,range,len);
    z = linspace(-range,range,len);
    [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(X,Y,Z);
    th = pi/2 - th;

    %% Calculation
    switch wavetype
    case 'x plane polarized'
        Ex = E0*exp(1i*cos(th).*(r-Rs)*k);
        Jx = Ex;
    case 'left circular'
     disp('not implemented yet');
    case 'right circular'
     disp('not implemented yet');
    case 'beam'
     disp('not implemented yet');
    otherwise
     disp(strcat(wavetype, ' is not a valid wave type'));
    end
    
end


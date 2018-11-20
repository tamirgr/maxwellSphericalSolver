function plotYlm()
    l = 5;
    m=1;
    len = 60;
        
    range =2;
    %x = linspace(-range,range,len);
    %y = linspace(-range,range,len);
    %z = linspace(-range,range,len);
    figure;
    [x,y,z] = sphere(len);
 %   [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(x,y,z);
    th = pi/2 -th;
    ylm = Ylm(l,m,th,phi);
    surf(x,y,z,real(ylm));
        axis equal
        shading interp
        colormap jet
end
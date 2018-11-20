function plotXlm()
    l = 2;
    m=1;
    len = 60;
        
    range =2;
    %x = linspace(-range,range,len);
    %y = linspace(-range,range,len);
    %z = linspace(-range,range,len);
    
    [x,y,z] = sphere(len);
 %   [X,Y,Z] = meshgrid(x,y,z);
    [phi,th,r] = cart2sph(x,y,z);
    th = pi/2 -th;
    [xlmTh,xlmPhi] = Xlm(th,phi,l,m);
    figure(1);
    
    surf(x,y,z,real(xlmTh));
    title('Theta');
    axis equal
    shading interp
    colormap jet
    figure(2);
    
    surf(x,y,z,real(xlmPhi));
    title('Phi');
    axis equal
        shading interp
        colormap jet
end
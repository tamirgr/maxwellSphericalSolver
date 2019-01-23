function displayFields( F1 , F2 , F3 ,X,Y,Z, n,l,m,disp)
 

if max(F1)
    figure;
    colormap('jet');
    slice(X,Y,Z,F1,disp,disp,disp);
    colorbar();
    shading interp
    %  caxis([-1    1]*1*10^-12);
    title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));
end

if max(F2)
figure;
colormap('jet');
slice(X,Y,Z,F2,disp,disp,disp);
colorbar();
shading interp
title(sprintf('EyReal n=%d, l=%d, m=%d',n,l,m));
end

if max(F3)
figure;
colormap('jet');
slice(X,Y,Z,F3,disp,disp,disp);
colorbar();
shading interp
title(sprintf('EzReal n=%d, l=%d, m=%d',n,l,m));
end
end
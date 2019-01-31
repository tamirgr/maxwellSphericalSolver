function displayFields( F1 , F2 , F3 ,X,Y,Z, n,l,m,disp)
 
axval = 6;

if max(F1)
    figure;
    colormap('jet');
    slice(X,Y,Z,F1,disp,disp,disp);
    colorbar();
    shading interp
    caxis([-axval    axval]);
    title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));
end

if max(F2)
figure;
colormap('jet');
slice(X,Y,Z,F2,disp,disp,disp);
colorbar();
shading interp
    caxis([-axval    axval]);
title(sprintf('EyReal n=%d, l=%d, m=%d',n,l,m));
end

if max(F3)
figure;
colormap('jet');
slice(X,Y,Z,F3,disp,disp,disp);
colorbar();
shading interp
    caxis([-axval    axval]);
title(sprintf('EzReal n=%d, l=%d, m=%d',n,l,m));
end
end
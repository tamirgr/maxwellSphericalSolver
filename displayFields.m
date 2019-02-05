function displayFields( F1 , F2 , F3 ,X,Y,Z, n,l,m,disp, level)
 
axval = 6;

if max(F1)
    figure;
    colormap('jet');
    slice(X,Y,Z,F1,disp,disp,disp);
    colorbar();
    shading interp
    caxis([-axval    axval]);
    title(sprintf('ExReal n=%d, l=%d, m=%d',n,l,m));
    switch level
        case 1
        movegui('northwest')
        case 2
        movegui('west')
        case 3
        movegui('southwest')
    end
end

if max(F2)
    figure;
    colormap('jet');
    slice(X,Y,Z,F2,disp,disp,disp);
    colorbar();
    shading interp
        caxis([-axval    axval]);
    title(sprintf('EyReal n=%d, l=%d, m=%d',n,l,m));
    switch level
        case 1
        movegui('north')
        case 2
        movegui('center')
        case 3
        movegui('south')
    end
end

if max(F3)
figure;
colormap('jet');
slice(X,Y,Z,F3,disp,disp,disp);
colorbar();
shading interp
    caxis([-axval    axval]);
title(sprintf('EzReal n=%d, l=%d, m=%d',n,l,m));
    switch level
        case 1
        movegui('northeast')
        case 2
        movegui('east')
        case 3
        movegui('southeast')
    end
end
end
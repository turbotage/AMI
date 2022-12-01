function line_filtered = full_line_filter(tvi)
    xmax = size(tvi,1);
    ymax = size(tvi,2);
    zmax = size(tvi,3);
    line_filtered=zeros(xmax,ymax,zmax);
    for i=1:zmax
        tvitemp = mat2gray(tvi(:,:,i));
        tvi1 = imbinarize(tvitemp,'adaptive','Sensitivity',0.1);
        tvi2 = line_masker(tvi1, true);
        tvi1 = line_masker(tvi1, false);
        tvi2 = line_masker2(tvi2, tvi1);
        line_filtered(:,:,i) = regionfill(tvi(:,:,i),tvi2);
    end
end
function varsum = var_map(dat, n_variation_chunks)
    zmax = size(dat,3);
    steplen = round(zmax/n_variation_chunks);
    under_step = 1;
    varsum = zeros(size(dat,1),size(dat,2));
    for i=steplen:steplen:zmax
        varsum = varsum + std(dat(:,:,under_step:i), 0, 3);
        under_step = under_step + steplen;
    end
end
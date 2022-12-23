function draw_std(mat, perc)
    if nargin < 2
        perc = [1,99];
    end
    std_mat = std(mat,[],3);
    perc = prctile(std_mat,perc);

    figure;
    imagesc(std_mat);
    clim(perc);
end
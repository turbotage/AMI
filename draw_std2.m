function draw_std2(mat1,mat2,perc)
    if nargin < 3
        perc = [1,99];
    end

    std_mat1 = std(mat1,[],3);
    std_mat2 = std(mat2,[],3);
    perc1 = prctile(std_mat1,perc,[1,2]);
    perc2 = prctile(std_mat2,perc,[1,2]);

    figure;
    subplot(1,2,1);
    imagesc(std_mat1);
    clim(perc1);

    subplot(1,2,2);
    imagesc(std_mat2);
    clim(perc2);
end
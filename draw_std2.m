function draw_std2(mat1,mat2,perc)
    if nargin < 3
        perc = [1,99];
    end

    std_mat1 = std(mat1,[],3);
    std_mat2 = std(mat2,[],3);
    perc1 = prctile(std_mat1,perc,[1,2]);
    perc2 = prctile(std_mat2,perc,[1,2]);

    t=tiledlayout(1,2),'TileSpacing','Compact');
	nexttile;
    imagesc(std_mat1);
    title('Low pass')
    clim(perc1);

    nexttile;
    imagesc(std_mat2);
    title('High pass')
    clim(perc2);

    title(t,'Standard deviation mask of H-scan mean channel intensity')
	cb = colorbar;
	cb.Layout.Tile = 'east';
end
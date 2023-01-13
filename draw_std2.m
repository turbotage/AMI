function draw_std2(chosen_dataset,mat1,mat2,perc,roimat)
% Rebeccas adaption of draw_std2 for H-scan. 
% Input:    chosen_dataset  String specifying the selected dataset for the 
%                           heading of its plot
%           mat1,mat2       r and b-encoded color channels
%           perc            resolution (has a default)
%           roimat          the binary mask specifying the roi

if nargin < 3
    perc = [1,99];
end

std_mat1 = std(mat1,[],3);
std_mat2 = std(mat2,[],3);
perc1 = prctile(std_mat1,perc,[1,2]);
perc2 = prctile(std_mat2,perc,[1,2]);

t=tiledlayout(1,3,'TileSpacing','Compact');
nexttile;
imagesc(roimat);
title('ROI')

nexttile;
imagesc(std_mat1);
title('Low pass')
clim(perc1);

nexttile;
imagesc(std_mat2);
title('High pass')
clim(perc2);

title(t,['Standard deviation mask of H-scan mean channel\newlineintensity for ' ...
    ,chosen_dataset])
cb = colorbar;
cb.Layout.Tile = 'east';
end
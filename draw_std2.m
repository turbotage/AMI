function draw_std2(chosen_dataset,mat1,mat2,perc,roimat)
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


% Get coordinates of the rois bdrys
% bdry = bwboundaries(roimat);
% bCoords=bdry{1};
% x = bCoords(:, 2);
% y = bCoords(:, 1);
% 
% topLine = min(x);
% bottomLine = max(x);
% leftColumn = min(y);
% rightColumn = max(y);
% width = bottomLine - topLine + 1;
% height = rightColumn - leftColumn + 1;
% 
% % std_mat1 cropped roi
% MaskedImg1 = std_mat1;
% MaskedImg1(~roimat) = 0;
% 
% croppedImg1 = imcrop(MaskedImg1, [topLine, leftColumn, width, height]);
% 
% % Display cropped image.
% nexttile;
% imagesc(croppedImg1);
% title('Low pass ROI');
% 
% % std_mat2 cropped roi
% MaskedImg2 = std_mat2;
% MaskedImg2(~roimat) = 0;
% 
% croppedImg2 = imcrop(MaskedImg2, [topLine, leftColumn, width, height]);
% 
% % Display cropped image.
% nexttile;
% imagesc(croppedImg2);
% title('High pass ROI');

title(t,['Standard deviation mask of H-scan mean channel\newlineintensity for ' ...
    ,chosen_dataset])
cb = colorbar;
cb.Layout.Tile = 'east';
end
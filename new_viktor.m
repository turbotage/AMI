%% load TVI
tvi_dat = matfile('181023_1311tvi_rsf.mat').tvi_dat;
%tvi_dat = single(matfile('181023_1311tvi_rs.mat').tvi_downsampled);
%% load SVD
svd_dat = matfile('181023_1311_rsf.mat').svd_dat;
svd_dat = sqrt(abs(hilbert(squeeze(svd_dat))));
%%
tvi_dat_f = mybandpass(tvi_dat,4,[0.5,8],500);
tvi_dat_f = full_line_filter(tvi_dat_f);
tvi_dat_f = medfilt3(tvi_dat_f);

%%
draw_pic(tvi_dat_f);

%%
figure;
std_tvi_img = std(tvi_dat_f, 0, 3);

imagesc(std_tvi_img);
colormap gray;

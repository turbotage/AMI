%% load TVI
h = matfile('181023_1311tvi_rs.mat');
tvi_dat = single(h.tvi_downsampled(:,:,1:1000));
%tvi_dat = sqrt(abs(hilbert(squeeze(tvi_dat(:,:,:)))));
clear h;

%% load SVD
svd_dat = single(load('181023_1311_rs.mat').rfmat_downsampled);
svd_dat = sqrt(abs(hilbert(squeeze(svd_dat(:,:,:)))));

%% SVD
%% Filter sections for SVD
svd_dat_filtered = svd_dat;
%% Frequency peak removal
svd_dat_filtered = svd_filter(svd_dat_filtered, true, false, false, false, true);
%% Butter
svd_dat_filtered = svd_filter(svd_dat_filtered, false, true, false, false, false);
%% Hard line removal
svd_dat_filtered = svd_filter(svd_dat_filtered, false, false, true, false, false);
%% General whole volume filter
svd_dat_filtered = svd_filter(svd_dat_filtered, false, false, false, true, false);
%% Per temporal slice convolution
svd_dat_filtered = svd_filter(svd_dat_filtered, false, false, false, false, true);
%%
draw_pic2(svd_dat_filtered, svd_dat);
%%
svd_hr = run_svd(svd_dat_filtered,0.002,0.1);

%%
draw_pic(svd_hr);

%%
svdsum = var_map(svd_hr,4);
figure;
imagesc(svdsum);
%axis('square');
%clim([0,max(svdsum,[],[1,2])]);

%% TVI
%% Filter sections for TVI
tvi_dat_filtered = tvi_dat;
%% Frequency peak removal
tvi_dat_filtered = tvi_filter(tvi_dat_filtered, true, false, false, false, true);
%% Butter
tvi_dat_filtered = tvi_filter(tvi_dat_filtered, false, true, false, false, false);
%% Hard line removal
tvi_dat_filtered = tvi_filter(tvi_dat_filtered, false, false, true, false, false);
%% General whole volume filter
tvi_dat_filtered = tvi_filter(tvi_dat_filtered, false, false, false, true, false);
%% Per temporal slice convolution
tvi_dat_filtered = tvi_filter(tvi_dat_filtered, false, false, false, false, true);

%%
draw_pic2(tvi_dat, tvi_dat_filtered);
%%
tvisum = var_map(tvi_dat_filtered,4);
figure;
imagesc(tvisum);

%axis('square');
%clim([0,0.01]);



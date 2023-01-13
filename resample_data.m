%% load TVI
tvi_dat = matfile('181023_0913tvi.mat').TVI_MAT;
disp('Ran load TVI')
%% load SVD
svd_dat = matfile('181023_0913.mat').rf_data_set;
disp('Ran load SVD')
%% load ROI
roi = matfile('181023_0913_ROI.mat').ROI;
%%
tvi_dat = tvi_dat(1:1000,:,:);
tvi_dat = single(myresample(tvi_dat));
save('181023_0913tvi_rsf.mat','tvi_dat');
%%
svd_dat = svd_dat(1:1000,:,:);
rfdat = single(myresample(svd_dat));
save('181023_0913_rsf.mat','rfdat');
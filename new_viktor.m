%% load TVI
tvi_dat = matfile('181023_1311tvi_rsf.mat').tvi_dat;
%tvi_dat = single(matfile('181023_1311tvi_rs.mat').tvi_downsampled);
%% load SVD
svd_dat = matfile('181023_1311_rsf.mat').svd_dat;
%%
tgc_vec = linspace(1,10,size(svd_dat,1))';
tgc_mat = repmat(tgc_vec,[1,size(svd_dat,2)]);
svd_dat = svd_dat .* tgc_mat;

%%
% tgc_vec = linspace(1,10,size(tvi_dat,1))';
% tgc_mat = repmat(tgc_vec,[1,size(tvi_dat,2)]);
% tvi_dat = tvi_dat .* tgc_mat;
% clear tgc_vec;
% clear tgc_mat;

%%
% figure;
% imagesc(svd_dat(:,:,10));
% clim([-6300, 6300]);
% 
% figure;
% imagesc(svd_dat_tgc(:,:,10));
% clim([-6300, 6300]);
%%
for i=1:size(svd_dat,3)
    svd_dat(:,:,i) = sqrt(abs(hilbert(squeeze(svd_dat(:,:,i)))));
end
%%
tvi_dat_f = mybandpass(tvi_dat,4,[3,30],500);
tvi_dat_f = full_line_filter(tvi_dat_f);
%%

%%
svd_dat_f = mybandpass(svd_dat,4,[3, 30],500);

%%
svd_dat_f = mymedfilt(svd_dat_f, [20,3]);

%%
svd_dat_fsvd = run_svd(svd_dat_f, 0.004, 0.90);
%%
std_svd_img = std(svd_dat_fsvd, 0, 3);
ptile = prctile(std_svd_img,[1,99]);

figure;
imagesc(std_svd_img(1:1000,:));
colormap gray;
clim(ptile);

%%

%%
draw_pic(svd_dat_fsvd);

%%
draw_pic(tvi_dat_f);

%%
figure;
std_tvi_img = std(tvi_dat_f, 0, 3);

imagesc(std_tvi_img);
colormap gray;

%%
%svd_conved = mypertempconv(svd_dat_fsvd, ones(20,3)/30);
svd_conved = mymedfilt(svd_dat_fsvd,[20,3]);

%%
svd_conved = mypertempconv(svd_dat_f, ones(20,3)/30);

%%
std_svd_img = std(svd_conved, 0, 3);
ptile = prctile(std_svd_img,[1,99]);

figure;
imagesc(std_svd_img(1:1000,:));
colormap gray;
clim(ptile);


%%

upper_signal_tvi = squeeze(mean(tvi_dat_f(300:305,60:65,:), [1,2]));
lower_signal_tvi = squeeze(mean(tvi_dat_f(800:805,60:65,:), [1,2]));

figure;
plot(upper_signal_tvi, 'r-');
hold on;
plot(lower_signal_tvi, 'b-');

%%
upper_signal_svd = squeeze(mean(svd_dat_fsvd(300:305,60:65,:), [1,2]));
lower_signal_svd = squeeze(mean(svd_dat_fsvd(800:805,60:65,:), [1,2]));

figure;
plot(upper_signal_svd, 'r-');
hold on;
plot(lower_signal_svd, 'b-');

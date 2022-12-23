%% load TVI
tvi_dat = matfile('181023_1311tvi_rsf.mat').tvi_dat;
disp('Ran load TVI')
%% load SVD
svd_dat = matfile('181023_1311_rsf.mat').rfdat;
disp('Ran load SVD')
%% TGC
tgc_vec = linspace(1,5,size(svd_dat,1))';
tgc_mat = repmat(tgc_vec,[1,size(svd_dat,2)]);
svd_dat = svd_dat .* tgc_mat;
clear tgc_mat;
clear tgc_vec;
disp('Applied TGC')
%% Hilbert
for i=1:size(svd_dat,3)
    svd_dat(:,:,i) = sqrt(abs(hilbert(squeeze(svd_dat(:,:,i)))));
end
clear i;
disp('Applied Hilbert')
%% Detrend
tvi_dat = tvi_dat - mean(tvi_dat,[1,2,3]);
svd_dat = svd_dat - mean(svd_dat,[1,2,3]);
disp('Ran Detrend')
%% TVI Filter
tvi_dat_f = mybandpass(tvi_dat,4,[3,30],500);
tvi_dat_f = mymedfilt(tvi_dat_f, [10,3]);
%tvi_dat_f = full_line_filter(tvi_dat_f);
disp('TVI Filter')
%%

%% Pre SVD Filters
svd_dat_f = mybandpass(svd_dat,4,[3, 30],500);
svd_dat_f = mymedfilt(svd_dat_f, [10,3]);
disp('Pre SVD Filters')
%% Run SVD
svd_dat_fsvd = run_svd(svd_dat_f, 5, 200);
disp('Run SVD')
%%
draw_std2(svd_dat_f, svd_dat_fsvd);
%%
draw_pic2(svd_dat_f, svd_dat_fsvd);

%%
draw_pic(tvi_dat_f);

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

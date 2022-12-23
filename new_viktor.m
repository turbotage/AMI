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


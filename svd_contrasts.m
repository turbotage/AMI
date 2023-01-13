%% load SVD
svd_dat = matfile('181023_1555_rsf.mat').rfdat;
disp('Ran load SVD')
% bandpass [3,20] for 1311 and 1555, [0.5,10] for 0913
bp_interval = [3,20];
%% load ROI
roi = matfile('181023_1555_ROI.mat').ROI;
disp('Ran load roi')
%% Hilbert
for i=1:size(svd_dat,3)
    svd_dat(:,:,i) = sqrt(abs(hilbert(squeeze(svd_dat(:,:,i)))));
end
clear i;
disp('Applied Hilbert')
%% Pre SVD Filters
svd_dat_f = mybandpass(svd_dat,4,bp_interval,500);
svd_dat_f = mymedfilt(svd_dat_f, [10,3]);
disp('Pre SVD Filters')
%%
optimal_interval_start = 0.75;
optimal_interval_len = 50;

start = optimal_interval_start * size(svd_dat_f,3);
svd_dat_fsvd = run_svd(svd_dat_f, start, start + optimal_interval_len);
disp('Run SVD')
%%
contrast = roi_metric(svd_dat_fsvd,roi);
%%
svd_dat_fsvd = mybandpass(svd_dat_fsvd,4,bp_interval,500);
%%
contrast_bp = roi_metric(svd_dat_fsvd,roi);


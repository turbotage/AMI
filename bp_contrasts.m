%% load SVD
svd_dat = matfile('181023_0913_rsf.mat').rfdat;
disp('Ran load SVD')
% bandpass [3,20] for 1311 and 1555, [0.5,10] for 0913
bp_interval = [0.5,10];
%% load ROI
roi = matfile('181023_0913_ROI.mat').ROI;
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
contrast = roi_metric(svd_dat_f,roi);
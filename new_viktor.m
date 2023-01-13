%% load TVI
tvi_dat = matfile('181023_1555tvi_rsf.mat').tvi_dat;
disp('Ran load TVI')
%% load SVD
svd_dat = matfile('181023_1555_rsf.mat').rfdat;
disp('Ran load SVD')
%% load ROI
roi = matfile('181023_1555_ROI.mat').ROI;
disp('Ran load roi')
%% TGC
tgc_vec = linspace(1,3,size(svd_dat,1))';
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
%% TVI Filter
% bandpass [3,20] for 1311 and 1555, [0.5,10] for 0913
tvi_dat_f = mybandpass(tvi_dat,4,[3,20],500);
tvi_dat_f = mymedfilt(tvi_dat_f, [10,3]);
%tvi_dat_f = full_line_filter(tvi_dat_f);
disp('TVI Filter')
%% Pre SVD Filters
% bandpass [3,20] for 1311 and 1555, [0.5,10] for 0913
svd_dat_f = mybandpass(svd_dat,4,[3,20],500);
svd_dat_f = mymedfilt(svd_dat_f, [10,3]);
disp('Pre SVD Filters')
%% Sweep SVD

singular_starts = [0,0.001,0.005,0.01,0.02,0.04,0.07,0.1,0.15,0.20,0.25,0.30,0.35,0.40,...
    0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95];
inters = [0,1,2,5,10,20,50,100,300,500];

ratios = svd_sweep(svd_dat_f, roi, singular_starts, inters);

disp('Ran sweep SVD')
%% Run SVD

start = 0.75 * size(svd_dat_f,3);
svd_dat_fsvd = run_svd(svd_dat_f, start, start + 50);
disp('Run SVD')

%% create TVI std maps and SVD std maps

svd_std = std(svd_dat_fsvd,[],3);
tvi_std = std(tvi_dat_f,[],3);

perc_svd = prctile(svd_std,[1,99],[1,2]);
perc_tvi = prctile(tvi_std,[1,99],[1,2]);

figure;
subplot(1,3,1);
imagesc(roi);
title('ROI');

subplot(1,3,2);
imagesc(svd_std);
clim(perc_svd);
title('SVD');

subplot(1,3,3);
imagesc(tvi_std);
clim(perc_tvi);
title('TVI');

